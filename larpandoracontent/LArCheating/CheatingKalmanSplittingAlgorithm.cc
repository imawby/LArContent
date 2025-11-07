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
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanFit.h"
#include "larpandoracontent/LArUtility/KalmanFilter.h"

#include "larpandoracontent/LArCheating/CheatingKalmanSplittingAlgorithm.h"

#include <numeric>

using namespace pandora;

namespace lar_content
{

CheatingKalmanSplittingAlgorithm::MCContaminant::MCContaminant(const pandora::CartesianVector &contStartPosition, 
    const pandora::CartesianVector &contEndPosition, const pandora::CartesianVector &startPosition, const pandora::CartesianVector &endPosition, 
    const pandora::CartesianVector &startDirection, const pandora::CartesianVector &endDirection) :
    m_contStartPosition(contStartPosition),
    m_contEndPosition(contEndPosition),
    m_startPosition(startPosition),
    m_endPosition(endPosition),
    m_startDirection(startDirection),
    m_endDirection(endDirection)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

CheatingKalmanSplittingAlgorithm::CheatingKalmanSplittingAlgorithm() :
    m_mcParticleListName("Input"), 
    m_secVertexListName("SecondaryVertices3D"),
    m_minClusterHits(50),
    m_minTargetMCHits(5),
    m_minFractionMerged(0.5f),
    m_slidingWindow(20),

    m_lBinSize(0.5f),
    m_endpointBuffer(2.f),
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
    std::map<const Cluster *, int> contaminantCounts;
    ClusterToSplitPositionsMap clusterToSplitPositionsMap;
    this->FillPandoraMaps(pClusterList, pCaloHitList, clusterToSplitPositionsMap, contaminantCounts);

    this->ProbeContaminants(pClusterList, pCaloHitList, pSecVertexList, clusterToSplitPositionsMap, contaminantCounts);


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingKalmanSplittingAlgorithm::FillPandoraMaps(const ClusterList *const pClusterList, const CaloHitList *const pCaloHitList, 
    ClusterToSplitPositionsMap &clusterToSplitPositionsMap, std::map<const Cluster *, int> &contaminantCounts)
{
    HitToMCParticleMap hitToMainMCParticleMap; // Dominating MCParticle
    std::map<const CaloHit*, std::vector<const MCParticle*>> hitToContMCParticleMap; // Contributing MCParticle

    // Fill our hit map
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        // Identify contributing MCParticles
        const MCParticleWeightMap hitMCParticleWeightMap(pCaloHit->GetMCParticleWeightMap());
        for (const auto &entry : hitMCParticleWeightMap)
        {
            if (entry.second > 0.1)
                hitToContMCParticleMap[pCaloHit].push_back(entry.first);
        }
        // Identify dominant MCParticles
        try
        {
            const MCParticle *const pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            hitToMainMCParticleMap[pCaloHit] = pMCParticle;
        }
        catch (...) { continue; }
    }

    // Now understand each cluster composition
    for (const Cluster *const pCluster : *pClusterList)
    {
        CaloHitList clusterHits;
        LArClusterHelper::GetAllHits(pCluster, clusterHits);

        // Fill contaminant hit lists
        MCParticleToHitListMap mainHitListMap, contHitListMap;
        for (const CaloHit *const pCaloHit : clusterHits)
        {
            if (hitToContMCParticleMap.find(pCaloHit) != hitToContMCParticleMap.end())
                for (const MCParticle *const pContMCParticle : hitToContMCParticleMap.at(pCaloHit))
                    contHitListMap[pContMCParticle].push_back(pCaloHit);

            if (hitToMainMCParticleMap.find(pCaloHit) != hitToMainMCParticleMap.end())
                mainHitListMap[hitToMainMCParticleMap.at(pCaloHit)].push_back(pCaloHit);
        }

        // Record number of contaminants (for tree filling later)
        int nContaminants(0);
        for (auto &entry : mainHitListMap)
            if (entry.second.size() >= 3)
                ++nContaminants;
        contaminantCounts.insert(std::make_pair(pCluster, nContaminants));

        // Find contaminants
        ContaminantMap contaminantMap;
        this->FindContaminants(mainHitListMap, contHitListMap, contaminantMap);

        // Now find split positions
        CartesianPointVector splitPositions;
        this->FindSplitPositions(pCluster, contaminantMap, splitPositions);

        // Add to map
        clusterToSplitPositionsMap.insert(std::make_pair(pCluster, splitPositions));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingKalmanSplittingAlgorithm::FindContaminants(const MCParticleToHitListMap &mainHitListMap, 
    const MCParticleToHitListMap &contHitListMap, ContaminantMap &contaminantMap)
{
    // Sort
    MCParticleVector mcParticleVector;
    for (const auto &entry : mainHitListMap) mcParticleVector.push_back(entry.first);
    std::sort(mcParticleVector.begin(), mcParticleVector.end(), PointerLessThan<MCParticle>());

    // Find contaminants
    for (const MCParticle *const pMCContaminant : mcParticleVector)
    {
        // Get hits MCParticle contributed to
        if (contHitListMap.find(pMCContaminant) == contHitListMap.end())
            continue;

        const CaloHitList &mainHits(mainHitListMap.at(pMCContaminant));
        const CaloHitList &contHits(contHitListMap.at(pMCContaminant));
        const HitType hitType(mainHits.front()->GetHitType());

        // Can we see an isolated particle? 
        unsigned int hitCount(0);
        for (const CaloHit *const pCaloHit : mainHits)
        {
            const MCParticleWeightMap hitMCParticleWeightMap(pCaloHit->GetMCParticleWeightMap());
            for (const auto &entry : hitMCParticleWeightMap)
            {
                if (entry.second > 0.9f)
                {
                    ++hitCount;
                    break;
                }
            }
        }
        if (hitCount < m_minTargetMCHits)
            continue;

        // Now build contaminant object
        const CartesianVector trueStart(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCContaminant->GetVertex(), hitType));
        // This isn't great for photons, but I think that's okay, they're not the target here
        const CartesianVector trueEnd(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCContaminant->GetEndpoint(), hitType));

        CartesianVector contStartPosition(0.f,0.f,0.f), contEndPosition(0.f,0.f,0.f);
        CartesianVector startPosition(0.f,0.f,0.f), endPosition(0.f,0.f,0.f);
        CartesianVector startDirection(0.f,0.f,0.f), endDirection(0.f,0.f,0.f);
        CartesianPointVector fitPositions;

        float contStartSepSq(std::numeric_limits<float>::max()), contEndSepSq(std::numeric_limits<float>::max());
        float startSepSq(std::numeric_limits<float>::max()), endSepSq(std::numeric_limits<float>::max());
        
        for (const CaloHit *const pContHit : contHits)
        {
            float this_startSepSq((pContHit->GetPositionVector() - trueStart).GetMagnitudeSquared());
            float this_endSepSq((pContHit->GetPositionVector() - trueEnd).GetMagnitudeSquared());

            if (this_startSepSq < contStartSepSq)
            {
                contStartSepSq = this_startSepSq;
                contStartPosition = pContHit->GetPositionVector();
            }

            if (this_endSepSq < contEndSepSq)
            { 
                contEndSepSq = this_endSepSq;
                contEndPosition = pContHit->GetPositionVector();
            }

            // Do fit for contributing positions
            fitPositions.push_back(pContHit->GetPositionVector());

            // Now find start/end positions used for splitting
            if (std::find(mainHits.begin(), mainHits.end(), pContHit) != mainHits.end())
            {
                if (this_startSepSq < startSepSq)
                {
                    startSepSq = this_startSepSq;
                    startPosition = pContHit->GetPositionVector();
                }

                if (this_endSepSq < endSepSq)
                { 
                    endSepSq = this_endSepSq;
                    endPosition = pContHit->GetPositionVector();
                }
            }
        }

        try
        {
            // Now direction.
            const TwoDSlidingFitResult clusterFit(&fitPositions, m_slidingWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));
            float startL(-1.f), startT(-1.f), endL(-1.f), endT(-1.f);
            clusterFit.GetLocalPosition(startPosition, startL, startT);
            clusterFit.GetLocalPosition(endPosition, endL, endT);
            clusterFit.GetGlobalFitDirection(startL, startDirection);
            clusterFit.GetGlobalFitDirection(endL, endDirection);

            // Now form and add to map
            contaminantMap.insert(std::make_pair(pMCContaminant, 
                MCContaminant(contStartPosition, contEndPosition, startPosition, endPosition, startDirection, endDirection)));
        }
        catch (...)
        {
            if (pMCContaminant->GetMomentum().GetMagnitude() < std::numeric_limits<float>::epsilon())
                continue;

            startDirection = LArGeometryHelper::ProjectDirection(this->GetPandora(), pMCContaminant->GetMomentum().GetUnitVector(), hitType);
            endDirection = startDirection;

            // Now form and add to map
            contaminantMap.insert(std::make_pair(pMCContaminant, 
                MCContaminant(contStartPosition, contEndPosition, startPosition, endPosition, startDirection, endDirection)));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingKalmanSplittingAlgorithm::FindSplitPositions(const Cluster *const pCluster, const ContaminantMap &contaminantMap, CartesianPointVector &splitPositions)
{
    // If there is only one MCParticle
    if (contaminantMap.size() < 2)
        return;
 
    MCParticleList contaminantList;
    for (auto &entry : contaminantMap)
        contaminantList.push_back(entry.first);

    try
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        const TwoDSlidingFitResult clusterFit(pCluster, m_slidingWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));
        const CartesianVector clusterMin(clusterFit.GetGlobalMinLayerPosition());
        const CartesianVector clusterMax(clusterFit.GetGlobalMaxLayerPosition());

        for (const MCParticle *const pThisContaminant : contaminantList)
        {
            if (contaminantMap.find(pThisContaminant) == contaminantMap.end())
                continue;

            const MCContaminant thisContaminant(contaminantMap.at(pThisContaminant));

            CartesianPointVector contPositions({thisContaminant.m_contStartPosition, thisContaminant.m_contEndPosition});
            CartesianPointVector positions({thisContaminant.m_startPosition, thisContaminant.m_endPosition});
            float l1(-1.f), t1(-1.f), l2(-1.f), t2(-1.f);
            clusterFit.GetLocalPosition(positions.at(0), l1, t1);
            clusterFit.GetLocalPosition(positions.at(1), l2, t2);
            CartesianPointVector directions({thisContaminant.m_startDirection, thisContaminant.m_endDirection});
            IntVector rejected;

            // Too close to cluster min/max?
            for (int i = 0; i < 2; ++i)
            {
                const float minSep((clusterMin - positions.at(i)).GetMagnitude());
                const float maxSep((clusterMax - positions.at(i)).GetMagnitude());

                if ((minSep < m_endpointBuffer) || (maxSep < m_endpointBuffer))
                    rejected.push_back(i);
            }

            if (rejected.size() == 2)
                continue;

        //////////////////////////////////////////
        // Draw contaminant
            //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &positions.at(0), "START", BLUE, 2);
            //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &positions.at(1), "END", BLUE, 2);
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        //////////////////////////////////////////

            // Does split pos live inside another particle?
            bool contained(false);

            // const float currentMinL(std::min(l1, l2));
            // const float currentMaxL(std::max(l1, l2));
            for (const MCParticle *const pTestContaminant : contaminantList)
            {
                if (pTestContaminant == pThisContaminant)
                    continue;

                if (contaminantMap.find(pTestContaminant) == contaminantMap.end())
                    continue;

                const MCContaminant testContaminant(contaminantMap.at(pTestContaminant));
                // ATTN: use contributing endpoints here!
                const CartesianVector &testStart(testContaminant.m_contStartPosition);
                const CartesianVector &testEnd(testContaminant.m_contEndPosition);

        //////////////////////////////////////////
        // Draw contaminant
                //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testStart, "START", RED, 2);
            //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testEnd, "END", RED, 2);
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        //////////////////////////////////////////
                float test_l1(-1.f), test_t1(-1.f), test_l2(-1.f), test_t2(-1.f);
                clusterFit.GetLocalPosition(testStart, test_l1, test_t1);
                clusterFit.GetLocalPosition(testEnd, test_l2, test_t2);
                const float testMinL(std::min(test_l1, test_l2));
                const float testMaxL(std::max(test_l1, test_l2));

                const bool startInside((l1 > testMinL) && (l1 < testMaxL));
                const bool endInside((l2 > testMinL) && (l2 < testMaxL));

                contained = (startInside && endInside);

                if (contained)
                    break;
            }

            if (contained)
                continue;

            // Reject collinear
            bool collinear(false);
            for (const MCParticle *const pTestContaminant : contaminantList)
            {
                if (pTestContaminant == pThisContaminant)
                    continue;

                if (contaminantMap.find(pTestContaminant) == contaminantMap.end())
                    continue;

                const MCContaminant testContaminant(contaminantMap.at(pTestContaminant));
                CartesianPointVector testPositions({testContaminant.m_startPosition, testContaminant.m_endPosition});
                CartesianPointVector testDirections({testContaminant.m_startDirection, testContaminant.m_endDirection});

                for (int iCurrent = 0; iCurrent < 2; ++iCurrent)
                {
                    if (std::find(rejected.begin(), rejected.end(), iCurrent) != rejected.end())
                        continue;

                    for (int iTest = 0; iTest < 2; ++iTest)
                    {
                        const float endpointSep((positions.at(iCurrent) - testPositions.at(iTest)).GetMagnitude());
                        float openingAngle(directions.at(iCurrent).GetOpeningAngle(testDirections.at(iTest)));
                        openingAngle *= (180.f / 3.14);

                        if ((endpointSep < 3.f) && ((openingAngle < 5.f) || (openingAngle > 175.f)))
                            collinear = true;

                        if (collinear)
                            break;
                    }

                    if (collinear)
                        rejected.push_back(iCurrent);
                }
            }

            // Finally, add in split positions
            for (int i = 0; i < 2; ++i)
            {
                if (std::find(rejected.begin(), rejected.end(), i) == rejected.end())
                    splitPositions.push_back(positions.at(i));
            }
        }
    }
    catch (...)
    {
        return;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingKalmanSplittingAlgorithm::ProbeContaminants(const ClusterList *const pClusterList, const CaloHitList *const pCaloHitList, 
    const VertexList *const pSecVertexList, const ClusterToSplitPositionsMap &clusterToSplitPositionsMap, 
    const std::map<const Cluster *, int> &contaminantCounts)
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

        // Does it have any contamination?
        const bool isContaminated((clusterToSplitPositionsMap.find(pCluster) != clusterToSplitPositionsMap.end()) && 
                                  (!clusterToSplitPositionsMap.at(pCluster).empty()));

        // Work out who truly owns this cluster
        int highestNHits(-1), matchedPDG(-1);
        std::map<const MCParticle*, CaloHitList> mcParticleHitListMap;

        for (const CaloHit *const pCaloHit : clusterHits)
        {
            try
            {
                const MCParticle *const pMainMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
                mcParticleHitListMap[pMainMCParticle].push_back(pCaloHit);

                if (static_cast<int>(mcParticleHitListMap.at(pMainMCParticle).size()) > highestNHits)
                {
                    highestNHits = mcParticleHitListMap.at(pMainMCParticle).size();
                    matchedPDG = pMainMCParticle->GetParticleId();
                }
            }
            catch (...) { continue; }
        }

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
                for (const CartesianVector &splitPosition : clusterToSplitPositionsMap.at(pCluster))
                {
                    float thisL(0.f), thisT(0.f);
                    clusterFit.GetLocalPosition(splitPosition, thisL, thisT);

                    vertexDrift.push_back(splitPosition.GetX());
                    vertexWire.push_back(splitPosition.GetZ());
                    vertexL.push_back(thisL);
                }
            }

            //////////////////////////////////
            // Pathway variables
            //////////////////////////////////
            HitKDTree2D kdTree_event;
            this->BuildKDTree(pCluster, pCaloHitList, kdTree_event);

            float cumulativeEnergy(0.f);

            const CaloHit *pPrevPathCaloHit(nullptr);

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
                gapSep.push_back(this->GetDistanceToGap(pPrevPathCaloHit, pPathCaloHit));
                secvertex.push_back(this->GetDistanceToSecVertex(pPathCaloHit, pSecVertexList, hitType));

                pPrevPathCaloHit = pPathCaloHit;
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

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "NContaminants", contaminantCounts.at(pCluster)));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ClusterPDG", std::abs(pCluster->GetParticleId())));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "BacktrackedPDG", matchedPDG));
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

float CheatingKalmanSplittingAlgorithm::GetDistanceToGap(const CaloHit *const pPrevHit, const CaloHit *const pCurrentHit) const
{
    if (!pPrevHit)
        return 1.f;

    const LArCaloHit *const pPrevLArHit(dynamic_cast<const LArCaloHit *>(pPrevHit));
    const LArCaloHit *const pCurrentLArHit(dynamic_cast<const LArCaloHit *>(pCurrentHit));

    if (!pPrevLArHit || !pCurrentLArHit)
        return 1.f;

    unsigned int prevTPCID(pPrevLArHit->GetLArTPCVolumeId());
    unsigned int currentTPCID(pCurrentLArHit->GetLArTPCVolumeId());

    if (prevTPCID != currentTPCID)
        return 0.f;

    unsigned int prevChildVolID(pPrevLArHit->GetDaughterVolumeId());
    unsigned int currentChildVolID(pCurrentLArHit->GetDaughterVolumeId());

    if (prevChildVolID != currentChildVolID)
        return 0.f;

    return 1.f;
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
