/**
 *  @file   larpandoracontent/LArCheating/CheatingKalmanSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the cheating cluster splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArKalmanHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanFit.h"
#include "larpandoracontent/LArUtility/KalmanFilter.h"

#include "larpandoracontent/LArCheating/CheatingKalmanSplittingAlgorithm.h"

#include <numeric>

using namespace pandora;

namespace lar_content
{

CheatingKalmanSplittingAlgorithm::CheatingKalmanSplittingAlgorithm() :
    m_mcParticleListName("Input"), 
    m_secVertexListName("SecondaryVertices3D"),
    m_minFractionMerged(0.5),
    m_minSecVertexAccuracy(5.f),
    m_writeFile(true),
    m_treeName("tree"),
    m_fileName("CheatingKalmanSplitting.root")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CheatingKalmanSplittingAlgorithm::~CheatingKalmanSplittingAlgorithm()
{
    if (m_writeFile)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingKalmanSplittingAlgorithm::Run()
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

    // Fill Pandora maps
    MCParticleToHitListMap mcParticleToHitListMap;
    HitToMCParticleMap hitToMCParticleMap;
    ClusterToMCParticleMap clusterToMCParticleMap;
    ClusterToMCParticleListMap clusterToMCParticleListMap;
    MCParticleSecVertexMap mcParticleSecVertexMap;

    this->FillPandoraMaps(pClusterList, pCaloHitList, pMCParticleList, pSecVertexList, mcParticleToHitListMap, hitToMCParticleMap, 
                          clusterToMCParticleMap, clusterToMCParticleListMap, mcParticleSecVertexMap);

    this->ProbeContaminants(pClusterList, pSecVertexList, clusterToMCParticleMap, clusterToMCParticleListMap, mcParticleToHitListMap);


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingKalmanSplittingAlgorithm::FillPandoraMaps(const ClusterList *const pClusterList, const CaloHitList *const pCaloHitList, const MCParticleList *const pMCParticleList, const VertexList *const pSecVertexList, MCParticleToHitListMap &mcParticleToHitListMap, HitToMCParticleMap &hitToMCParticleMap, ClusterToMCParticleMap &clusterToMCParticleMap, ClusterToMCParticleListMap &clusterToMCParticleListMap, MCParticleSecVertexMap &mcParticleSecVertexMap)
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

            if (totalMCHits < 5)
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

void CheatingKalmanSplittingAlgorithm::ProbeContaminants(const ClusterList *const pClusterList, const VertexList *const pSecVertexList, 
    ClusterToMCParticleMap &clusterToMCParticleMap, ClusterToMCParticleListMap &clusterToMCParticleListMap, MCParticleToHitListMap &mcParticleToHitListMap)
{
    ClusterList clusterList(*pClusterList);

    // For tree
    int clusterCount(0);

    for (const Cluster *const pCluster : clusterList)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if (clusterToMCParticleMap.find(pCluster) == clusterToMCParticleMap.end())
            continue;

        if (clusterToMCParticleListMap.find(pCluster) == clusterToMCParticleListMap.end())
            continue;

        const MCParticle *const pBestMatch(clusterToMCParticleMap.at(pCluster));

        KalmanFit kalmanFit(LArKalmanHelper::PerformKalmanFit(this->GetPandora(), pCluster, 5.f));

        if (kalmanFit.m_positions.empty())
            continue;

        // Does it have any contamination?
        if ((clusterToMCParticleListMap.at(pCluster).size() == 1) && (clusterToMCParticleListMap.at(pCluster).front() == pBestMatch))
            continue;

        // Make a fit for the cluster
        try
        {
            const TwoDSlidingFitResult clusterFit(&kalmanFit.m_positions, 20, LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));
            const CartesianVector clusterMin(clusterFit.GetGlobalMinLayerPosition());
            const CartesianVector clusterMax(clusterFit.GetGlobalMaxLayerPosition());

            std::vector<float> vertexDrift, vertexWire, vertexL;
            std::vector<float> driftCoord, wireCoord, lCoord, tCoord;
            std::vector<int> trackID, hitPDG, isInPath;
            std::vector<float> longitudinal, energy, angle, secvertex;

            CaloHitList clusterHits;
            LArClusterHelper::GetAllHits(pCluster, clusterHits);

            for (const CaloHit *const pCaloHit : clusterHits)
            {
                driftCoord.push_back(pCaloHit->GetPositionVector().GetX());
                wireCoord.push_back(pCaloHit->GetPositionVector().GetZ());
                float thisHitL(0.f), thisHitT(0.f);
                clusterFit.GetLocalPosition(pCaloHit->GetPositionVector(), thisHitL, thisHitT);
                lCoord.push_back(thisHitL);
                tCoord.push_back(thisHitT);
                isInPath.push_back((std::find(kalmanFit.m_caloHitList.begin(), kalmanFit.m_caloHitList.end(), pCaloHit) == 
                                    kalmanFit.m_caloHitList.end()) ? 0 : 1);

                try
                {
                    const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
                    hitPDG.push_back(pMCParticle->GetParticleId());
                    trackID.push_back((size_t)(intptr_t *)pMCParticle->GetUid());
                }
                catch (...)
                {
                    hitPDG.push_back(-1);
                    trackID.push_back(-1);
                }
            }

            //////////////////////////////////
            // Visualise contaminant
            ClusterList visualiseClusters({pCluster});
            PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &visualiseClusters, "Cluster", RED));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &clusterMin, "Fit", RED, 2));
            PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &clusterMax, "Fit", RED, 2));
            //////////////////////////////////

            for (const MCParticle *const pMCContaminant : clusterToMCParticleListMap.at(pCluster))
            {
                // if (pMCContaminant == pBestMatch)
                //     continue;

                if (mcParticleToHitListMap.find(pMCContaminant) == mcParticleToHitListMap.end())
                    throw;

                CartesianVector trueVertex(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCContaminant->GetVertex(), hitType));
                float thisVertexL(0.f), thisVertexT(0.f);
                clusterFit.GetLocalPosition(trueVertex, thisVertexL, thisVertexT);
                vertexDrift.push_back(trueVertex.GetX());
                vertexWire.push_back(trueVertex.GetZ());
                vertexL.push_back(thisVertexL);

                //////////////////////////////////
                CartesianVector trueJam(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCContaminant->GetVertex(), hitType));
                PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &trueJam, "True Vertex", BLACK, 2));
                PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
                //////////////////////////////////
            }

            ++clusterCount;

            //////////////////////////////////
            PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            //////////////////////////////////

            for (unsigned int i = 0; i < kalmanFit.m_positions.size(); ++i)
            {
                try 
                {
                    float rL(0.f), rT(0.f);
                    clusterFit.GetLocalPosition(kalmanFit.m_positions.at(i), rL, rT);
                    longitudinal.push_back(rL);
                }
                catch (...)
                {
                    continue;
                }

                energy.push_back(kalmanFit.m_caloHitList.at(i)->GetElectromagneticEnergy());

                if ((i == 0) || (i == 1))
                {
                    angle.push_back(-1.f);
                }
                else
                {
                    try
                    { 
                        angle.push_back(kalmanFit.m_directions.at(i).GetOpeningAngle(kalmanFit.m_directions.at(i-1)));
                    }
                    catch (...)
                    {
                        angle.push_back(-1.f);
                    }
                }

                float bestSep(-1.f);

                for (const Vertex *const pSecVertex : *pSecVertexList)
                {
                    const CartesianVector secVtxPos(LArGeometryHelper::ProjectPosition(this->GetPandora(), pSecVertex->GetPosition(), hitType));
                    const float thisSep((secVtxPos - kalmanFit.m_positions.at(i)).GetMagnitude());

                    if (bestSep < 0.f)
                    {
                        bestSep = thisSep;
                    }
                    else
                    {
                        bestSep = std::min(bestSep, thisSep);
                    }
                }

                secvertex.push_back(bestSep);
            }


            //////////////////////////////////
                ///            std::vector<float> driftCoord, wireCoord, lCoord, tCoord, trackID, hitPDG, isInPath;


            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ClusterCount", clusterCount));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "DriftCoord", &driftCoord));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "WireCoord", &wireCoord));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "LCoord", &lCoord));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "TCoord", &tCoord));

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "TrackID", &trackID));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "HitPDG", &hitPDG));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "IsInPath", &isInPath));

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "VertexDift", &vertexDrift));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "VertexWire", &vertexWire));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "VertexL", &vertexL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Energy", &energy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Longitudinal", &longitudinal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Angle", &angle));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecVertex", &secvertex));
            std::cout << "YO YO YO" << std::endl;
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
            //////////////////////////////////
        }
        catch(...)
        {
            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingKalmanSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "WriteFile", m_writeFile));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
