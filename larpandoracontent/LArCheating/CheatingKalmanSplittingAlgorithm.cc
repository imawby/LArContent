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
    //PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
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

    this->ProbeContaminants(pClusterList, pCaloHitList, pSecVertexList, clusterToMCParticleMap, clusterToMCParticleListMap, mcParticleToHitListMap);


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

void CheatingKalmanSplittingAlgorithm::FindPath(const Cluster *const pCluster, const TwoDSlidingFitResult &clusterFit, 
    std::map<int, std::pair<const CaloHit*, float>> &clusterPath)
{
    // Get cluster hits
    CaloHitList clusterHits;
    LArClusterHelper::GetAllHits(pCluster, clusterHits);

    // Fit l decomposition
    const float binSize(0.5f);
    std::map<int, std::vector<std::pair<const CaloHit*, float>>> lDecomposition;

    for (const CaloHit *const pCaloHit : clusterHits)
    {
        float thisHitL(0.f), thisHitT(0.f);
        clusterFit.GetLocalPosition(pCaloHit->GetPositionVector(), thisHitL, thisHitT);
        lDecomposition[std::floor(thisHitL / binSize)].push_back(std::make_pair(pCaloHit, thisHitT));
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

void CheatingKalmanSplittingAlgorithm::PerformKalmanFit(const std::map<int, std::pair<const CaloHit*, float>> &clusterPath, 
    std::vector<float> &posDiffDist, std::vector<float> &energyDiffDist, std::vector<float> &widthDiffDist, std::vector<float> &scatterDist, 
    std::vector<float> &mahalanobisDist)
{
    // Get total energy
    float totalEnergy(0.f);
    for (const auto &entry : clusterPath)
        totalEnergy += entry.second.first->GetElectromagneticEnergy();

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
    CartesianVector previousDirection(0.f, 0.f, 0.f);

    for (const auto &entry : clusterPath)
    {
        if (!skippedFirst)
        {
            skippedFirst = true;

            posDiffDist.push_back(-1.f);
            energyDiffDist.push_back(-1.f);
            widthDiffDist.push_back(-1.f);
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
        const float hitWidthDiff(-1.f);//std::fabs(pCaloHit->GetCellSize1() - tempState(3)));

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
        widthDiffDist.push_back(hitWidthDiff);
        scatterDist.push_back(openingAngleL);
        mahalanobisDist.push_back(mahalanobisDistance);

        // Catch direction 
        previousDirection = newDirection;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingKalmanSplittingAlgorithm::ProbeContaminants(const ClusterList *const pClusterList, const CaloHitList *const pCaloHitList, 
    const VertexList *const pSecVertexList, ClusterToMCParticleMap &clusterToMCParticleMap, ClusterToMCParticleListMap &clusterToMCParticleListMap, 
    MCParticleToHitListMap &mcParticleToHitListMap)
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

        // Does it have any contamination?
        const MCParticle *const pBestMatch(clusterToMCParticleMap.at(pCluster));

        if ((clusterToMCParticleListMap.at(pCluster).size() == 1) && (clusterToMCParticleListMap.at(pCluster).front() == pBestMatch))
            continue;

        try
        {
            // Make a fit for the cluster
            const TwoDSlidingFitResult clusterFit(pCluster, 20, LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));
            const CartesianVector clusterMin(clusterFit.GetGlobalMinLayerPosition());
            const CartesianVector clusterMax(clusterFit.GetGlobalMaxLayerPosition());

            // Find pathway through the cluster
            std::map<int, std::pair<const CaloHit*, float>> clusterPath;
            this->FindPath(pCluster, clusterFit, clusterPath);
                      
            if (clusterPath.empty())
                continue;

            // Get total energy
            float totalEnergy(0.f);
            for (const auto &entry : clusterPath)
                totalEnergy += entry.second.first->GetElectromagneticEnergy();

            // Now do Kalman fit
            std::vector<float> posDiffDist, energyDiffDist, widthDiffDist, scatterDist, mahalanobisDist;
            this->PerformKalmanFit(clusterPath, posDiffDist, energyDiffDist, widthDiffDist, scatterDist, mahalanobisDist);

            ////////////////////////////////////////////////////////
            // Fill out other tree stuff
            std::vector<float> vertexDrift, vertexWire, vertexL;
            //std::vector<float> driftCoord, wireCoord, lCoord, tCoord;
            //std::vector<int> trackID, hitPDG, isInPath;
            std::vector<float> longitudinal, transverse, energy, secvertex, width, gapSep, eventHitSep, clusterHitSep;

            // First plotting stuff
            // CaloHitList clusterHits;
            // LArClusterHelper::GetAllHits(pCluster, clusterHits);

            // for (const CaloHit *const pCaloHit : clusterHits)
            // {
            //     driftCoord.push_back(pCaloHit->GetPositionVector().GetX());
            //     wireCoord.push_back(pCaloHit->GetPositionVector().GetZ());
            //     float thisHitL(0.f), thisHitT(0.f);
            //     clusterFit.GetLocalPosition(pCaloHit->GetPositionVector(), thisHitL, thisHitT);
            //     lCoord.push_back(thisHitL);
            //     tCoord.push_back(thisHitT);

            //     try
            //     {
            //         const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            //         hitPDG.push_back(pMCParticle->GetParticleId());
            //         trackID.push_back((size_t)(intptr_t *)pMCParticle->GetUid());
            //     }
            //     catch (...)
            //     {
            //         hitPDG.push_back(-1);
            //         trackID.push_back(-1);
            //     }
            // }

            //////////////////////////////////
            // // Visualise contaminant
            // ClusterList visualiseClusters({pCluster});
            // PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &visualiseClusters, "Cluster", RED));
            // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &clusterMin, "Fit", RED, 2));
            // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &clusterMax, "Fit", RED, 2));
            //////////////////////////////////

            // Then contaminant vertices
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
                // CartesianVector trueJam(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCContaminant->GetVertex(), hitType));
                // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &trueJam, "True Vertex", BLACK, 2));
                // PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
                //////////////////////////////////
            }

            ++clusterCount;

            //////////////////////////////////
            //PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            //////////////////////////////////

            // Then pathway stuff
            //std::map<int, std::pair<const CaloHit*, float>> clusterPath;

            float cumEnergy(0.f);

            for (const auto &entry : clusterPath)
            {
                cumEnergy += (entry.second.first->GetElectromagneticEnergy() / totalEnergy);
                longitudinal.push_back(entry.first);
                transverse.push_back(entry.second.second);
                energy.push_back(cumEnergy);
                width.push_back(entry.second.first->GetCellSize1());
                gapSep.push_back(this->GetDistanceToGap(entry.second.first->GetPositionVector()));

                // Get separation within cluster
                float minClusterSepSq(std::numeric_limits<float>::max());
                bool skip(true);
                bool isClusterSepSet(false);
                for (const auto &entry2 : clusterPath)
                {
                    if (entry.second.first == entry2.second.first)
                    {
                        skip = false;
                        continue;
                    }

                    if (skip)
                        continue;

                    isClusterSepSet = true;
                    const float thisClusterSepSq((entry.second.first->GetPositionVector() - entry2.second.first->GetPositionVector()).GetMagnitudeSquared());
                    minClusterSepSq = std::min(thisClusterSepSq, minClusterSepSq);
                }

                clusterHitSep.push_back(isClusterSepSet ? std::sqrt(minClusterSepSq) : -1.f);

                // Get cluster hits
                CaloHitList clusterHits;
                LArClusterHelper::GetAllHits(pCluster, clusterHits);
                eventHitSep.push_back(this->GetDistanceToEventHit(entry.second.first, clusterHits, pCaloHitList));

                float bestSep(-1.f);

                for (const Vertex *const pSecVertex : *pSecVertexList)
                {
                    const CartesianVector secVtxPos(LArGeometryHelper::ProjectPosition(this->GetPandora(), pSecVertex->GetPosition(), hitType));
                    const float thisSep((secVtxPos - entry.second.first->GetPositionVector()).GetMagnitude());

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
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ClusterCount", clusterCount));
            // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "DriftCoord", &driftCoord));
            // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "WireCoord", &wireCoord));
            // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "LCoord", &lCoord));
            // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "TCoord", &tCoord));
            // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "TrackID", &trackID));
            // PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "HitPDG", &hitPDG));

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
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "WidthDiffDist", &widthDiffDist));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "MahalanobisDist", &mahalanobisDist));
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

float CheatingKalmanSplittingAlgorithm::GetDistanceToEventHit(const CaloHit *const pCaloHit, const CaloHitList &clusterHits, 
    const CaloHitList *const pEventHits)
{
    float minDistSq(std::numeric_limits<float>::max());

    for (const CaloHit *const pEventHit : *pEventHits)
    {
        if (std::find(clusterHits.begin(), clusterHits.end(), pEventHit) != clusterHits.end())
            continue;

        minDistSq = std::min(minDistSq, (pCaloHit->GetPositionVector() - pEventHit->GetPositionVector()).GetMagnitudeSquared());
    }

    return (pEventHits->empty() ? -1.f : std::sqrt(minDistSq));
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
