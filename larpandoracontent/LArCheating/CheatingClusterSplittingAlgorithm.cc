/**
 *  @file   larpandoracontent/LArCheating/CheatingClusterSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the cheating cluster splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArCheating/CheatingClusterSplittingAlgorithm.h"

#include <numeric>

using namespace pandora;

namespace lar_content
{

CheatingClusterSplittingAlgorithm::CheatingClusterSplittingAlgorithm() :
    m_mcParticleListName("Input"), 
    m_secVertexListName("SecondaryVertices3D"),
    m_minFractionMerged(0.5),
    m_minSecVertexAccuracy(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingClusterSplittingAlgorithm::Run()
{
    //////////////////////////////////
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    //////////////////////////////////

    // Get view hits
    const CaloHitList *pCaloHitList(nullptr);
    if (PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_SUCCESS;

    if ((!pCaloHitList) || pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    // Get view clusters
    const ClusterList *pClusterList(nullptr);
    if (PandoraContentApi::GetList(*this, m_clusterListName, pClusterList) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_SUCCESS;

    if ((!pClusterList) || pClusterList->empty())
        return STATUS_CODE_SUCCESS;

    // Get MCParticles
    const MCParticleList *pMCParticleList(nullptr);
    if (PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_SUCCESS;

    if ((!pMCParticleList) || pMCParticleList->empty())
        return STATUS_CODE_SUCCESS;

    // Get Secondary vertices
    const VertexList *pSecVertexList(nullptr);
    if (PandoraContentApi::GetList(*this, m_secVertexListName, pSecVertexList) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_SUCCESS;

    if ((!pSecVertexList) || pSecVertexList->empty())
        return STATUS_CODE_SUCCESS;

    bool changed(true);

    while (changed)
    {
        //std::cout << "new loop" << std::endl;
        changed = false;

        // Fill Pandora maps
        MCParticleToHitListMap mcParticleToHitListMap;
        HitToMCParticleMap hitToMCParticleMap;
        ClusterToMCParticleMap clusterToMCParticleMap;
        ClusterToMCParticleListMap clusterToMCParticleListMap;
        MCParticleSecVertexMap mcParticleSecVertexMap;

        this->FillPandoraMaps(pClusterList, pCaloHitList, pMCParticleList, pSecVertexList, mcParticleToHitListMap, hitToMCParticleMap, 
                              clusterToMCParticleMap, clusterToMCParticleListMap, mcParticleSecVertexMap);

        changed = this->SplitClusters(pClusterList, mcParticleSecVertexMap, clusterToMCParticleMap, clusterToMCParticleListMap, 
                                      mcParticleToHitListMap);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingClusterSplittingAlgorithm::FillPandoraMaps(const ClusterList *const pClusterList, const CaloHitList *const pCaloHitList, const MCParticleList *const pMCParticleList, const VertexList *const pSecVertexList, MCParticleToHitListMap &mcParticleToHitListMap, HitToMCParticleMap &hitToMCParticleMap, ClusterToMCParticleMap &clusterToMCParticleMap, ClusterToMCParticleListMap &clusterToMCParticleListMap, MCParticleSecVertexMap &mcParticleSecVertexMap)
{
    // Fill sec vertex map
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        const Vertex *pBestVertex(nullptr);
        float separationSq(std::numeric_limits<float>::max());

        for (const Vertex *const pSecVertex : *pSecVertexList)
        {
            const float thisSepSq((pSecVertex->GetPosition() - pMCParticle->GetVertex()).GetMagnitudeSquared());

            if (thisSepSq < separationSq)
            {
                separationSq = thisSepSq;
                pBestVertex = pSecVertex;
            }
        }

        if ((pBestVertex) && (std::sqrt(separationSq) < m_minSecVertexAccuracy))
            mcParticleSecVertexMap.insert(std::make_pair(pMCParticle, pBestVertex->GetPosition()));
    }

    // Fill our hit maps
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pMainMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            mcParticleToHitListMap[pMainMCParticle].push_back(pCaloHit);
            hitToMCParticleMap[pCaloHit] = pMainMCParticle;
        }
        catch (...) { continue; }
    }

    // Now understand each cluster composition
    for (const Cluster *const pCluster : *pClusterList)
    {
        std::unordered_map<const pandora::MCParticle *, FloatVector> clusterMCParticleToHitListMap;

        CaloHitList clusterHits;
        LArClusterHelper::GetAllHits(pCluster, clusterHits);

        for (const CaloHit *const pCaloHit : clusterHits)
        {
            if (hitToMCParticleMap.find(pCaloHit) == hitToMCParticleMap.end())
                continue;

            clusterMCParticleToHitListMap[hitToMCParticleMap.at(pCaloHit)].push_back(pCaloHit->GetElectromagneticEnergy());
        }

        // Find best match, and any match where more than 50% of the hits are in another particle
        int highestNHits(0);
        float highestEnergy(-1.f);
        const MCParticle *pBestMCParticle(nullptr);

        for (const auto &entry : clusterMCParticleToHitListMap)
        {
            const int nHits(entry.second.size());
            const float energySum(std::accumulate(entry.second.begin(), entry.second.end(), 0.f));

            if (nHits == highestNHits)
            {
                if (energySum > highestEnergy)
                {
                    highestNHits = entry.second.size();
                    highestEnergy = energySum;
                    pBestMCParticle = entry.first;
                }
            }
            else if (nHits > highestNHits)
            {
                highestNHits = entry.second.size();
                highestEnergy = energySum;
                pBestMCParticle = entry.first;
            }

            // is it more than 50%?
            // is it target reco? i.e. reconstructable?
            const int totalMCHits(mcParticleToHitListMap.find(entry.first) == mcParticleToHitListMap.end() ? 
                0 : mcParticleToHitListMap.at(entry.first).size());

            if (totalMCHits < 10)
                continue;

            const float particleCompleteness(totalMCHits == 0 ? 
                0.f : static_cast<float>(entry.second.size()) / static_cast<float>(totalMCHits));

            if (particleCompleteness > m_minFractionMerged)
                clusterToMCParticleListMap[pCluster].push_back(entry.first);
        }

        clusterToMCParticleMap[pCluster] = pBestMCParticle;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingClusterSplittingAlgorithm::SplitClusters(const ClusterList *const pClusterList, MCParticleSecVertexMap &mcParticleSecVertexMap, 
    ClusterToMCParticleMap &clusterToMCParticleMap, ClusterToMCParticleListMap &clusterToMCParticleListMap, MCParticleToHitListMap &mcParticleToHitListMap)
{
    bool madeChanges(false);
    ClusterList clusterList(*pClusterList);

    for (const Cluster *const pCluster : clusterList)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if (clusterToMCParticleMap.find(pCluster) == clusterToMCParticleMap.end())
            continue;

        if (clusterToMCParticleListMap.find(pCluster) == clusterToMCParticleListMap.end())
            continue;

        const MCParticle *const pBestMatch(clusterToMCParticleMap.at(pCluster));

        // Make a fit for the cluster
        try
        {
            const TwoDSlidingFitResult clusterFit(pCluster, 20, LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));
            const CartesianVector clusterMin(clusterFit.GetGlobalMinLayerPosition());
            const CartesianVector clusterMax(clusterFit.GetGlobalMaxLayerPosition());

            for (const MCParticle *const pMCContaminant : clusterToMCParticleListMap.at(pCluster))
            {
                if (pMCContaminant == pBestMatch)
                    continue;

                if (mcParticleToHitListMap.find(pMCContaminant) == mcParticleToHitListMap.end())
                    throw;

                if (mcParticleSecVertexMap.find(pMCContaminant) == mcParticleSecVertexMap.end())
                    continue;

                CartesianVector trueVertex(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCContaminant->GetVertex(), hitType));
                CartesianVector secVertex(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcParticleSecVertexMap.at(pMCContaminant), hitType));

                // we should make sure that the sec vertex isn't near the start of the cluster (because then we'd have nothing to split)
                const float minSep((clusterMin - secVertex).GetMagnitude());
                const float maxSep((clusterMax - secVertex).GetMagnitude());
                if ((minSep < m_minSecVertexAccuracy) || (maxSep < m_minSecVertexAccuracy))
                    continue;

                //////////////////////////////////
                // Visualise before split
                // ClusterList visualiseClusters({pCluster});
                // PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &visualiseClusters, "Cluster", RED));
                // CartesianVector jam(LArGeometryHelper::ProjectPosition(this->GetPandora(), mcParticleSecVertexMap.at(pMCContaminant), hitType));
                // CartesianVector trueJam(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCContaminant->GetVertex(), hitType));
                // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &jam, "SecVertex", BLUE, 2));
                // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &trueJam, "True Vertex", BLACK, 2));
                // PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
                //////////////////////////////////

                //std::cout << "Lets split the particle" << std::endl;
                //std::cout << "AAAAA" << std::endl;
                // First divide hits
                CaloHitList caloHitList1, caloHitList2;
                if (this->DivideCaloHits(clusterFit, trueVertex, caloHitList1, caloHitList2) != STATUS_CODE_SUCCESS)
                    continue;
                //std::cout << "BBBB" << std::endl;
                // Now we can do the splitting
                ClusterVector clusterSplittingList;
                if (this->SplitCluster(pCluster, caloHitList1, caloHitList2, clusterSplittingList) != STATUS_CODE_SUCCESS)
                    continue;
                //std::cout << "CCCCC" << std::endl;
                //////////////////////////////////
                // // Visualise after split
                // ClusterList visualiseClusters1({clusterSplittingList.at(0)});
                // ClusterList visualiseClusters2({clusterSplittingList.at(1)});
                // PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &visualiseClusters1, "Cluster1", GREEN));
                // PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &visualiseClusters2, "Cluster2", VIOLET));
                // float rL(0.f), rT(0.f);
                // clusterFit.GetLocalPosition(trueVertex, rL, rT);
                // CartesianVector visualiseSplitPoint(0.f, 0.f, 0.f);
                // clusterFit.GetGlobalPosition(rL, 0.f, visualiseSplitPoint);
                // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &visualiseSplitPoint, "SecVertex", BLUE, 2));
                // PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
                //////////////////////////////////
                //std::cout << "DDDD" << std::endl;
                madeChanges = true;
                break;
            }
        }
        catch(...)
        {
            continue;
        }
    }

    return madeChanges;
}

//------------------------------------------------------------------------------------------------------------------------------------------

// Stolen from TwoDSlidingFitSplittingAlgorithm (this alg should really inherit from that alg)
StatusCode CheatingClusterSplittingAlgorithm::DivideCaloHits(const TwoDSlidingFitResult &slidingFitResult, const CartesianVector &splitPosition, 
    CaloHitList &firstCaloHitList, CaloHitList &secondCaloHitList) const
{
    float rL(0.f), rT(0.f);
    slidingFitResult.GetLocalPosition(splitPosition, rL, rT);

    const Cluster *const pCluster(slidingFitResult.GetCluster());
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CaloHit *const pCaloHit = *hitIter;

            float thisL(0.f), thisT(0.f);
            slidingFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), thisL, thisT);

            if (thisL < rL)
            {
                firstCaloHitList.push_back(pCaloHit);
            }
            else
            {
                secondCaloHitList.push_back(pCaloHit);
            }
        }
    }

    if (firstCaloHitList.empty() || secondCaloHitList.empty())
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

// Stolen from ClusterSplittingAlgorithm, which TwoDSlidingFitSplittingAlgorithm inherits from
StatusCode CheatingClusterSplittingAlgorithm::SplitCluster(const Cluster *const pCluster, const CaloHitList &caloHitList1, const CaloHitList &caloHitList2,
    ClusterVector &clusterSplittingList) const
{
    // Split cluster into two CaloHit lists
    PandoraContentApi::Cluster::Parameters parameters1, parameters2;
    parameters1.m_caloHitList = caloHitList1;
    parameters2.m_caloHitList = caloHitList2;

    if (parameters1.m_caloHitList.empty() || parameters2.m_caloHitList.empty())
        return STATUS_CODE_NOT_ALLOWED;

    // Begin cluster fragmentation operations
    const ClusterList clusterList(1, pCluster);
    std::string clusterListToSaveName, clusterListToDeleteName;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                             PandoraContentApi::InitializeFragmentation(*this, clusterList, clusterListToDeleteName, clusterListToSaveName));

    // Create new clusters
    const Cluster *pCluster1(NULL), *pCluster2(NULL);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters1, pCluster1));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters2, pCluster2));

    clusterSplittingList.push_back(pCluster1);
    clusterSplittingList.push_back(pCluster2);

    // End cluster fragmentation operations
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingClusterSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "SecVertexListName", m_secVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MinFractionMerged", m_minFractionMerged));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MinSecVertexAccuracy", m_minSecVertexAccuracy));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
