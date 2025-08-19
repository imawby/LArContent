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
    m_minClusterHits(50),
    m_minTargetMCHits(5),
    m_minFractionMerged(0.5f),
    m_slidingWindow(20),
    m_lBinSize(0.5f),
    m_endpointBuffer(3.f),
    m_searchRegion1D(20.f),
    m_writeVisInfo(false),
    m_treeName("tree"),
    m_fileName("CheatingKalmanSplitting.root")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CheatingKalmanSplittingAlgorithm::~CheatingKalmanSplittingAlgorithm()
{
    PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingKalmanSplittingAlgorithm::Run()
{
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

    this->FillPandoraMaps(pClusterList, pCaloHitList, pMCParticleList, pSecVertexList, mcParticleToHitListMap, hitToMCParticleMap, 
                          clusterToMCParticleMap, clusterToMCParticleListMap);

    this->ProbeContaminants(pClusterList, pCaloHitList, pSecVertexList, clusterToMCParticleMap, clusterToMCParticleListMap, mcParticleToHitListMap);


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingKalmanSplittingAlgorithm::FillPandoraMaps(const ClusterList *const pClusterList, const CaloHitList *const pCaloHitList, const MCParticleList *const pMCParticleList, const VertexList *const pSecVertexList, MCParticleToHitListMap &mcParticleToHitListMap, HitToMCParticleMap &hitToMCParticleMap, ClusterToMCParticleMap &clusterToMCParticleMap, ClusterToMCParticleListMap &clusterToMCParticleListMap)
{
    // Fill our hit map
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
        std::unordered_map<const pandora::MCParticle *, CaloHitList> clusterMCParticleToHitListMap;

        CaloHitList clusterHits;
        LArClusterHelper::GetAllHits(pCluster, clusterHits);

        for (const CaloHit *const pCaloHit : clusterHits)
        {
            if (hitToMCParticleMap.find(pCaloHit) == hitToMCParticleMap.end())
                continue;

            clusterMCParticleToHitListMap[hitToMCParticleMap.at(pCaloHit)].push_back(pCaloHit);
        }

        // Find best match, and any MCParticle that have > 5 hits and 50% of them are in this cluster
        int highestNHits(0);
        float highestEnergy(-1.f);
        const MCParticle *pBestMCParticle(nullptr);

        for (const auto &entry : clusterMCParticleToHitListMap)
        {
            const int nHits(entry.second.size());
            float energySum(0.f);
            for (const CaloHit *const pCaloHit : entry.second)
                energySum += pCaloHit->GetElectromagneticEnergy();

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

            // Is there a significant contamination?
            if (nHits < m_minTargetMCHits)
                continue;

            // Is a lot of the MC hits missing?
            // const int totalMCHits(mcParticleToHitListMap.find(entry.first) == mcParticleToHitListMap.end() ? 
            //     0 : mcParticleToHitListMap.at(entry.first).size());

            // const float particleCompleteness(totalMCHits == 0 ? 
            //     0.f : static_cast<float>(entry.second.size()) / static_cast<float>(totalMCHits));

            // if (particleCompleteness < m_minFractionMerged)
            //     continue;

            clusterToMCParticleListMap[pCluster].push_back(entry);
        }

        clusterToMCParticleMap[pCluster] = pBestMCParticle;
    }
}


//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingKalmanSplittingAlgorithm::FindPath(const Cluster *const pCluster, const TwoDSlidingFitResult &clusterFit, 
    ClusterPath &clusterPath)
{
    // Get cluster hits
    CaloHitList clusterHits;
    LArClusterHelper::GetAllHits(pCluster, clusterHits);

    // Fit l decomposition
    std::map<int, std::vector<std::pair<const CaloHit*, float>>> lDecomposition;

    for (const CaloHit *const pCaloHit : clusterHits)
    {
        float thisHitL(0.f), thisHitT(0.f);
        clusterFit.GetLocalPosition(pCaloHit->GetPositionVector(), thisHitL, thisHitT);
        lDecomposition[std::floor(thisHitL / m_lBinSize)].push_back(std::make_pair(pCaloHit, thisHitT));
    }

    // Find path
    for (std::pair<int, std::vector<std::pair<const CaloHit*, float>>> entry : lDecomposition)
    {
        if (entry.second.size() == 1)
        {
            clusterPath.insert(std::make_pair(entry.first, entry.second.front()));
        }
        else
        {
            // pick pair with smallest t
            std::pair<const CaloHit*, float> pBestHit(entry.second.front());
            float smallestT(std::numeric_limits<float>::max());

            for (std::pair<const CaloHit*, float> ambiguousHit : entry.second)
            {
                if (ambiguousHit.second < smallestT)
                {
                    smallestT = ambiguousHit.second;
                    pBestHit = ambiguousHit;
                }
            }

            clusterPath.insert(std::make_pair(entry.first, pBestHit));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingKalmanSplittingAlgorithm::ProbeContaminants(const ClusterList *const pClusterList, const CaloHitList *const pCaloHitList, 
    const VertexList *const pSecVertexList, ClusterToMCParticleMap &clusterToMCParticleMap, ClusterToMCParticleListMap &clusterToMCParticleListMap, 
    MCParticleToHitListMap &mcParticleToHitListMap)
{
    ClusterList clusterList(*pClusterList);

    for (const Cluster *const pCluster : clusterList)
    {
        // Enough hits?
        CaloHitList clusterHits;
        LArClusterHelper::GetAllHits(pCluster, clusterHits);

        if (clusterHits.size() < m_minClusterHits)
            continue;

        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));

        if (clusterToMCParticleMap.find(pCluster) == clusterToMCParticleMap.end())
            continue;

        if (clusterToMCParticleListMap.find(pCluster) == clusterToMCParticleListMap.end())
            continue;

        // Does it have any contamination?
        const MCParticle *const pBestMatch(clusterToMCParticleMap.at(pCluster));

        // if ((clusterToMCParticleListMap.at(pCluster).size() == 1) && (clusterToMCParticleListMap.at(pCluster).front() == pBestMatch))
        //     continue;

        bool isContaminated(!((clusterToMCParticleListMap.at(pCluster).size() == 1) && (clusterToMCParticleListMap.at(pCluster).front().first == pBestMatch)));

        try
        {
            // Make a fit for the cluster
            const TwoDSlidingFitResult clusterFit(pCluster, m_slidingWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));
            const CartesianVector clusterMin(clusterFit.GetGlobalMinLayerPosition());
            const CartesianVector clusterMax(clusterFit.GetGlobalMaxLayerPosition());

            // Find pathway through the cluster
            ClusterPath clusterPath;
            this->FindPath(pCluster, clusterFit, clusterPath);
                      
            if (clusterPath.empty())
                continue;

            // Get total energy of path
            float totalEnergy(0.f);
            for (const auto &entry : clusterPath)
                totalEnergy += entry.second.first->GetElectromagneticEnergy();

            //////////////////////////////////
            // Variables to fill
            //////////////////////////////////
            std::vector<float> longitudinal, transverse, energy, secvertex, width, gapSep, eventHitSep, clusterHitSep;
            std::vector<float> posDiffDist, energyDiffDist, scatterDist, mahalanobisDist;
            std::vector<float> vertexDrift, vertexWire, vertexL;
            std::vector<float> driftCoord, wireCoord, lCoord, tCoord;
            std::vector<int> trackID, hitPDG, isInPath;

            //////////////////////////////////
            // Kalman fit variables
            //////////////////////////////////
            this->PerformKalmanFit(clusterPath, totalEnergy, posDiffDist, energyDiffDist, scatterDist, mahalanobisDist);

            //////////////////////////////////
            // Plotting variables
            //////////////////////////////////
            if (m_writeVisInfo)
            {
                for (const CaloHit *const pCaloHit : clusterHits)
                {
                    driftCoord.push_back(pCaloHit->GetPositionVector().GetX());
                    wireCoord.push_back(pCaloHit->GetPositionVector().GetZ());
                    float thisHitL(0.f), thisHitT(0.f);
                    clusterFit.GetLocalPosition(pCaloHit->GetPositionVector(), thisHitL, thisHitT);
                    lCoord.push_back(thisHitL);
                    tCoord.push_back(thisHitT);

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
            }

            //////////////////////////////////
            // Splitting positions
            //////////////////////////////////
            if (isContaminated)
            {
                for (const auto &entry : clusterToMCParticleListMap.at(pCluster))
                {
                    const MCParticle *const pMCContaminant(entry.first);
                    CartesianVector trueVertex(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCContaminant->GetVertex(), hitType));
                    const float minVertexSep((clusterMin - trueVertex).GetMagnitude());
                    const float maxVertexSep((clusterMax - trueVertex).GetMagnitude());

                    if ((minVertexSep < m_endpointBuffer) || (maxVertexSep < m_endpointBuffer))
                        continue;

                    // Find the closest hit to true position
                    const CaloHitList &contaminantHits(entry.second);
                    float smallestSep(std::numeric_limits<float>::max());
                    const CaloHit *pClosestHit(nullptr);

                    for (const CaloHit *const pContaminantHit : contaminantHits)
                    {
                        const float thisSepSq((pContaminantHit->GetPositionVector() - trueVertex).GetMagnitudeSquared());

                        if (thisSepSq < smallestSep)
                        {
                            smallestSep = thisSepSq;
                            pClosestHit = pContaminantHit;
                        }
                    }

                    if (!pClosestHit)
                        continue;

                    const float minSep((clusterMin - pClosestHit->GetPositionVector()).GetMagnitude());
                    const float maxSep((clusterMax - pClosestHit->GetPositionVector()).GetMagnitude());

                    if ((minSep < m_endpointBuffer) || (maxSep < m_endpointBuffer))
                        continue;

                    float thisVertexL(0.f), thisVertexT(0.f);
                    clusterFit.GetLocalPosition(pClosestHit->GetPositionVector(), thisVertexL, thisVertexT);

                    vertexDrift.push_back(trueVertex.GetX());
                    vertexWire.push_back(trueVertex.GetZ());
                    vertexL.push_back(thisVertexL);
                }
            }

            // Reset isContaminated i.e. do we have any split positions?
            isContaminated = !vertexL.empty();

            //////////////////////////////////
            // Pathway variables
            //////////////////////////////////
            HitKDTree2D kdTree_event;
            this->BuildKDTree(pCluster, pCaloHitList, kdTree_event);

            float cumulativeEnergy(0.f);

            for (ClusterPath::iterator iter = clusterPath.begin(); iter != clusterPath.end(); ++iter)
            {
                longitudinal.push_back(iter->first);
                transverse.push_back(iter->second.second);

                const CaloHit *const pPathCaloHit(iter->second.first);

                cumulativeEnergy += (pPathCaloHit->GetElectromagneticEnergy() / totalEnergy);
                energy.push_back(cumulativeEnergy);
                width.push_back(pPathCaloHit->GetCellSize1());
                eventHitSep.push_back(this->GetDistanceToEventHit(kdTree_event, pCluster, pPathCaloHit));
                clusterHitSep.push_back(this->GetDistanceToClusterHit(iter, clusterPath.end()));
                gapSep.push_back(this->GetDistanceToGap(pPathCaloHit->GetPositionVector()));
                secvertex.push_back(this->GetDistanceToSecVertex(pPathCaloHit, pSecVertexList, hitType));
            }

            //////////////////////////////////
            // Fill Tree
            //////////////////////////////////
            if (m_writeVisInfo)
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "DriftCoord", &driftCoord));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "WireCoord", &wireCoord));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "LCoord", &lCoord));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "TCoord", &tCoord));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "TrackID", &trackID));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "HitPDG", &hitPDG));
            }

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "IsContaminated", (isContaminated ? 1 : 0)));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "VertexDift", &vertexDrift));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "VertexWire", &vertexWire));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "VertexL", &vertexL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Energy", &energy));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "HitWidth", &width));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "GapSep", &gapSep));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "EventHitSep", &eventHitSep));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ClusterHitSep", &clusterHitSep));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Longitudinal", &longitudinal));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Transverse", &transverse));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "Angle", &scatterDist));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecVertex", &secvertex));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PosDiffDist", &posDiffDist));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "EnergyDiffDist", &energyDiffDist));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "MahalanobisDist", &mahalanobisDist));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
        }
        catch(...)
        {
            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingKalmanSplittingAlgorithm::PerformKalmanFit(const ClusterPath &clusterPath, 
    const float totalEnergy, std::vector<float> &posDiffDist, std::vector<float> &energyDiffDist, std::vector<float> &scatterDist, 
    std::vector<float> &mahalanobisDist)
{
    // Kalman Config
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const HitType view(clusterPath.begin()->second.first->GetHitType());
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
    const float m_kalmanDelta(1.f), m_kalmanProcessVarCoeff(1.f), m_kalmanMeasurementVarCoeff(1.f);
    const float processVariance{m_kalmanProcessVarCoeff * pitch * pitch};
    const float measurementVariance{m_kalmanMeasurementVarCoeff * pitch * pitch};
   
    // Initialise Kalman fit
    Eigen::VectorXd init(3);
    const CaloHit *const seedHit(clusterPath.begin()->second.first);
    const float seedL(clusterPath.begin()->first), seedT(clusterPath.begin()->second.second);
    init << seedL, seedT, (seedHit->GetElectromagneticEnergy() / totalEnergy);

    KalmanFilter3D kalmanFilter3D(m_kalmanDelta, processVariance, measurementVariance, init);

    bool skippedFirst(false);

    for (const auto &entry : clusterPath)
    {
        if (!skippedFirst)
        {
            skippedFirst = true;

            posDiffDist.push_back(-1.f);
            energyDiffDist.push_back(-1.f);
            scatterDist.push_back(-4.f);
            mahalanobisDist.push_back(-1.f);

            continue;
        }

        // Make next kalman step
        kalmanFilter3D.Predict();
        const KalmanFilter3D::StateVector &tempState(kalmanFilter3D.GetTemporaryState());

        // Compare
        const CartesianVector thisPosition(entry.first, 0.f, entry.second.second);
        const CartesianVector predPosition(tempState(0), 0.f, tempState(1));
        const float separation((thisPosition - predPosition).GetMagnitude());
        const CaloHit *const pCaloHit(entry.second.first);
        const float energyDiff(std::fabs((pCaloHit->GetElectromagneticEnergy() / totalEnergy) - tempState(2))); 

        // Update filter
        Eigen::VectorXd eigenXd(3);
        eigenXd << thisPosition.GetX(), thisPosition.GetZ(), (pCaloHit->GetElectromagneticEnergy() / totalEnergy);
        float mahalanobisDistance(kalmanFilter3D.GetMahalanobisDistance(eigenXd, true));
        kalmanFilter3D.Update(eigenXd);

        // Get scatter angle
        const CartesianVector newDirection(kalmanFilter3D.GetState()(3), 0.f, kalmanFilter3D.GetState()(4));
        float openingAngleL(-999.f);

        try
        {
            const float openingAngleT = CartesianVector(0.f, 0.f, 1.f).GetOpeningAngle(newDirection);
            openingAngleL = CartesianVector(1.f, 0.f, 0.f).GetOpeningAngle(newDirection);
            openingAngleL *= (openingAngleT > (M_PI * 0.5f)) ? (-1.f) : 1.f;
        }
        catch (...) {};

        // Store info
        posDiffDist.push_back(separation);
        energyDiffDist.push_back(energyDiff);
        scatterDist.push_back(openingAngleL);
        mahalanobisDist.push_back(mahalanobisDistance);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingKalmanSplittingAlgorithm::BuildKDTree(const Cluster *const pCluster, const CaloHitList *const pCaloHitList, 
    HitKDTree2D &kdTree_event)
{
    // KD tree for cluster...
    CaloHitList clusterHits;
    LArClusterHelper::GetAllHits(pCluster, clusterHits);

    CaloHitList eventHits(*pCaloHitList);

    for (const CaloHit *const pClusterHit : clusterHits)
        eventHits.remove(pClusterHit);

    HitKDNode2DList kdNode2DList_event;
    KDTreeBox kdTreeBox_event(fill_and_bound_2d_kd_tree(eventHits, kdNode2DList_event));

    kdTree_event.build(kdNode2DList_event, kdTreeBox_event);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CheatingKalmanSplittingAlgorithm::GetDistanceToEventHit(HitKDTree2D &kdTree_event, const Cluster *const pCluster, 
    const CaloHit *const pCaloHit)
{
    HitKDNode2DList found;
    KDTreeBox searchRegionHits(build_2d_kd_search_region(pCaloHit, m_searchRegion1D, m_searchRegion1D));
    kdTree_event.search(searchRegionHits, found);

    float minDistSq(std::numeric_limits<float>::max());
    for (const auto &hit : found)
        minDistSq = std::min(minDistSq, (pCaloHit->GetPositionVector() - hit.data->GetPositionVector()).GetMagnitudeSquared());

    return (found.empty() ? -1.f : std::sqrt(minDistSq));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CheatingKalmanSplittingAlgorithm::GetDistanceToClusterHit(const ClusterPath::iterator &currentHit, const ClusterPath::iterator &endIter)
{
    const ClusterPath::iterator nextHit(std::next(currentHit));

    if (nextHit == endIter)
        return -1.f;

    return (currentHit->second.first->GetPositionVector() - nextHit->second.first->GetPositionVector()).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CheatingKalmanSplittingAlgorithm::GetDistanceToGap(const CartesianVector &position2D) const
{
    const DetectorGapList detectorGapList(this->GetPandora().GetGeometry()->GetDetectorGapList());

    float minDist(std::numeric_limits<float>::max());

    if (detectorGapList.empty())
        minDist = -1.f;

    for (const DetectorGap *const pDetectorGap : detectorGapList)
    {
        const LineGap *const pLineGap(dynamic_cast<const LineGap *>(pDetectorGap));

        if (!pLineGap)
            continue;
        
        const LineGapType lineGapType(pLineGap->GetLineGapType());
            
        if (lineGapType == TPC_DRIFT_GAP)
        {
            minDist = std::min(std::fabs(pLineGap->GetLineStartX() - position2D.GetX()), minDist);
            minDist = std::min(std::fabs(pLineGap->GetLineEndX() - position2D.GetX()), minDist);
        }

        if ((lineGapType == TPC_WIRE_GAP_VIEW_U) || (lineGapType == TPC_WIRE_GAP_VIEW_V) || (lineGapType == TPC_WIRE_GAP_VIEW_W))
        {
            minDist = std::min(std::fabs(pLineGap->GetLineStartZ() - position2D.GetZ()), minDist);
            minDist = std::min(std::fabs(pLineGap->GetLineEndZ() - position2D.GetZ()), minDist);
        }
    }

    return minDist;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CheatingKalmanSplittingAlgorithm::GetDistanceToSecVertex(const CaloHit *const pCaloHit, const VertexList *const pSecVertexList, 
    const HitType hitType)
{
    float bestSepSq(std::numeric_limits<float>::max());

    for (const Vertex *const pSecVertex : *pSecVertexList)
    {
        const CartesianVector secVtxPos(LArGeometryHelper::ProjectPosition(this->GetPandora(), pSecVertex->GetPosition(), hitType));
        bestSepSq = std::min(bestSepSq, (secVtxPos - pCaloHit->GetPositionVector()).GetMagnitudeSquared());
    }

    return (pSecVertexList->empty() ? -1.f : std::sqrt(bestSepSq));
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
        XmlHelper::ReadValue(xmlHandle, "MinClusterHits", m_minClusterHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MinTargetMCHits", m_minTargetMCHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MinFractionMerged", m_minFractionMerged));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "SlidingWindow", m_slidingWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "LBinSize", m_lBinSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "EndpointBuffer", m_endpointBuffer));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "SearchRegion1D", m_searchRegion1D));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "WriteVisualisationInfo", m_writeVisInfo));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
