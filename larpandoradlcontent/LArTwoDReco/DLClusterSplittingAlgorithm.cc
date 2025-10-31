/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/DLClusterSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the two dimensional sliding fit splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArUtility/KalmanFilter.h"
#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"


#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArTwoDReco/DLClusterSplittingAlgorithm.h"

#include <torch/script.h>
#include <torch/torch.h>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

DLClusterSplittingAlgorithm::DLClusterSplittingAlgorithm() :
    m_secVertexListName("SecondaryVertices3D"),
    m_pSecVertexList(nullptr),
    m_minClusterHits(50),
    m_slidingWindow(20),
    m_lBinSize(0.5f),
    m_searchRegion1D(20.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLClusterSplittingAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

    // Get view hits
    const CaloHitList *pCaloHitList(nullptr);
    if (PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_SUCCESS;

    if ((!pCaloHitList) || pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    // Get secondary vertices (it's okay if the list is empty)
    PandoraContentApi::GetList(*this, m_secVertexListName, m_pSecVertexList);

    // Get clusters
    const ClusterList *pClusterList(nullptr);
    PandoraContentApi::GetList(*this, m_clusterListName, pClusterList);

    if ((!pClusterList) || pClusterList->empty())
        return STATUS_CODE_SUCCESS;

    // Make current...
    if (PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListName) != STATUS_CODE_SUCCESS)
    {
        std::cout << "isobel cannot make current - sad" << std::endl;
        throw;
    }

    ClusterVector internalClusterVector(pClusterList->begin(), pClusterList->end());
    //internalClusterVector.sort(LArClusterHelper::SortByNHits);

    // Probe clusters
    for (const Cluster *const pCluster : internalClusterVector)
        this->ThisDivideCaloHits(pCluster);

    m_pSecVertexList = nullptr;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLClusterSplittingAlgorithm::ThisDivideCaloHits(const Cluster *const pCluster)
{
    // Enough hits?
    CaloHitList clusterHits;
    LArClusterHelper::GetAllHits(pCluster, clusterHits);

    if (clusterHits.size() < m_minClusterHits)
        return STATUS_CODE_NOT_FOUND;

    // Can we make fit?
    try
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        const TwoDSlidingFitResult clusterFit(pCluster, m_slidingWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));

        // Find pathway through the cluster
        ClusterPath clusterPath;
        this->FindPath(pCluster, clusterFit, clusterPath);
                      
        if (clusterPath.empty())
            return STATUS_CODE_NOT_FOUND;

        // Get features
        Features features;
        this->FillFeatures(clusterPath, clusterFit, features);

        std::vector<int> clusterLPositions;
        for (auto entry1 : clusterPath)
            clusterLPositions.push_back(entry1.first);

        // Get windows
        IntVector splitIndices;
        this->GetWindows(features, splitIndices);

        // Remove any that are too close to the nu vertex (network picks up on hit sharing)
        FloatVector lSplit;

        //NeutrinoVertices3D
        const VertexList *pVertexList(nullptr);
        PandoraContentApi::GetList(*this, "NeutrinoVertices3D", pVertexList);

        if (!pVertexList || pVertexList->empty())
            return STATUS_CODE_SUCCESS;

        const CartesianVector nuVertexPosition(pVertexList->front()->GetPosition());

        //std::cout << "-----------------------------------------" << std::endl;
        //std::cout << "-----------------------------------------" << std::endl;
        for (const int splitIndex : splitIndices)
        {
            int clusterSplitIndex(clusterLPositions.at(splitIndex));
            CartesianVector position(clusterPath.at(clusterSplitIndex).first->GetPositionVector());

            if ((nuVertexPosition - position).GetMagnitude() > 5.f)
            {
                float thisL(0.f), thisT(0.f);
                clusterFit.GetLocalPosition(position, thisL, thisT);
                lSplit.push_back(thisL);
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "SplitPoint", VIOLET, 2);
            }
        }

        int nSplitPositions(lSplit.size());

        if (nSplitPositions == 0)
        {
            //PandoraMonitoringApi::ViewEvent(this->GetPandora());
            return STATUS_CODE_NOT_FOUND;
        }
    
        ClusterList visCluster({pCluster});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &visCluster, "Cluster", BLACK);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());

        // insert max numbers..
        lSplit.insert(lSplit.begin(), std::numeric_limits<float>::min());
        lSplit.insert(lSplit.end(), std::numeric_limits<float>::max());
        std::vector<CaloHitList> splitClusterHits((nSplitPositions + 1), CaloHitList());

        for (const CaloHit *const pCaloHit : clusterHits)
        {
            float thisL(0.f), thisT(0.f);
            clusterFit.GetLocalPosition(pCaloHit->GetPositionVector(), thisL, thisT);

            for (int i = 0; i <= nSplitPositions; ++i)
            {
                if ((thisL > lSplit.at(i)) && (thisL < lSplit.at(i + 1)))
                {
                    splitClusterHits[i].push_back(pCaloHit);
                    break;
                }
            }
        }

        // Split clusters
        // Begin cluster fragmentation operations
        const ClusterList clusterList(1, pCluster);
        std::string clusterListToSaveName, clusterListToDeleteName;

        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            PandoraContentApi::InitializeFragmentation(*this, clusterList, clusterListToDeleteName, clusterListToSaveName));

        for (CaloHitList &caloHitList : splitClusterHits)
        {
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList = caloHitList;

            const Cluster *pNewCluster(nullptr);
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewCluster));
        }

        // End cluster fragmentation operations
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));
    }
    catch (...)
    {
        return STATUS_CODE_NOT_FOUND;
    }

    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLClusterSplittingAlgorithm::FindPath(const Cluster *const pCluster, const TwoDSlidingFitResult &clusterFit, ClusterPath &clusterPath) const
{
    // Get cluster hits
    CaloHitList clusterHits;
    LArClusterHelper::GetAllHits(pCluster, clusterHits);

    // Fit l decomposition
    for (const CaloHit *const pCaloHit : clusterHits)
    {
        float thisHitL(0.f), thisHitT(0.f);
        clusterFit.GetLocalPosition(pCaloHit->GetPositionVector(), thisHitL, thisHitT);
        const int lBinIndex(std::floor(thisHitL / m_lBinSize));

        if (clusterPath.find(lBinIndex) != clusterPath.end())
            if (thisHitT > clusterPath.at(lBinIndex).second)
                continue;

        clusterPath[lBinIndex] = std::make_pair(pCaloHit, thisHitT);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLClusterSplittingAlgorithm::FillFeatures(const ClusterPath &clusterPath, const TwoDSlidingFitResult &clusterFit, Features &features) const
{
    // Get path hits, and their total energy
    float totalEnergy(0.f);
    CaloHitList clusterPathHits;
    for (const auto &entry : clusterPath)
    {
        const CaloHit *const pPathCaloHit(entry.second.first);
        totalEnergy += pPathCaloHit->GetElectromagneticEnergy();
        clusterPathHits.push_back(pPathCaloHit);
    }

    // Get KDTree
    const CaloHitList *pCaloHitList(nullptr);
    if (PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList) != STATUS_CODE_SUCCESS)
        return;

    if ((!pCaloHitList) || pCaloHitList->empty())
        return;

    HitKDTree2D kdTree;
    HitKDNode2DList kdNode2DList;
    KDTreeBox kdTreeBox(fill_and_bound_2d_kd_tree(*pCaloHitList, kdNode2DList));
    kdTree.build(kdNode2DList, kdTreeBox);

    // Kalman Config
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const HitType view(clusterPath.begin()->second.first->GetHitType());
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
    const float m_kalmanDelta(1.f), m_kalmanProcessVarCoeff(1.f), m_kalmanMeasurementVarCoeff(1.f);
    const float processVariance{m_kalmanProcessVarCoeff * pitch * pitch};
    const float measurementVariance{m_kalmanMeasurementVarCoeff * pitch * pitch};

    // Initialise Kalman fit
    Eigen::VectorXd init(2);
    const float seedL(clusterPath.begin()->first), seedT(clusterPath.begin()->second.second);
    init << seedL, seedT;
    KalmanFilter2D kalmanFilter2D(m_kalmanDelta, processVariance, measurementVariance, init);

    bool processedFirst(false);
    float cumulativeEnergy(0.f);

    // Fill features
    for (ClusterPath::const_iterator iter = clusterPath.begin(); iter != clusterPath.end(); ++iter)
    {
        if (!processedFirst)
        {
            features.m_theta.push_back(-4.f);
            processedFirst = true;
        }
        else
        {
            // Make next step
            kalmanFilter2D.Predict();
            // Update feature
            Eigen::VectorXd eigenXd(2);
            const CartesianVector thisPosition(iter->first, 0.f, iter->second.second);
            eigenXd << thisPosition.GetX(), thisPosition.GetZ();
            kalmanFilter2D.Update(eigenXd);
            // Get scatter angle
            const CartesianVector newDirection(kalmanFilter2D.GetState()(2), 0.f, kalmanFilter2D.GetState()(3));
            float openingAngleL(-4.f);
            try
            {
                const float openingAngleT = CartesianVector(0.f, 0.f, 1.f).GetOpeningAngle(newDirection);
                openingAngleL = CartesianVector(1.f, 0.f, 0.f).GetOpeningAngle(newDirection);
                openingAngleL *= (openingAngleT > (M_PI * 0.5f)) ? (-1.f) : 1.f;
            }
            catch (...) {};
            features.m_theta.push_back(openingAngleL);
        }

        const CaloHit *const pPathCaloHit(iter->second.first);
        cumulativeEnergy += pPathCaloHit->GetElectromagneticEnergy();

        features.m_transverse.push_back(iter->second.second);
        features.m_energy.push_back(cumulativeEnergy / totalEnergy);
        features.m_hitWidth.push_back(pPathCaloHit->GetCellSize1());
        features.m_secVertex.push_back(
            this->GetDistanceToSecVertex(pPathCaloHit, pPathCaloHit->GetHitType()));
        features.m_eventHitSep.push_back(
            this->GetDistanceToEventHit(kdTree, pPathCaloHit, clusterPathHits));
        features.m_clusterHitSep.push_back(
            this->GetDistanceToClusterHit(iter, clusterPath.end()));
        features.m_gapSep.push_back(
            this->GetDistanceToGap(pPathCaloHit->GetPositionVector(), pPathCaloHit->GetHitType()));
    }

    // Smooth
    features.SmoothFeature(features.m_transverse);
    features.SmoothFeature(features.m_energy);
    features.SmoothFeature(features.m_hitWidth);
    features.SmoothFeature(features.m_theta);
    features.SmoothFeature(features.m_secVertex);
    features.SmoothFeature(features.m_eventHitSep);
    features.SmoothFeature(features.m_clusterHitSep);
    features.SmoothFeature(features.m_gapSep);

    // Normalise -.-
    for (unsigned int i = 0; i < features.m_transverse.size(); ++i)
    {
        features.m_transverse[i] = features.NormaliseTransverse(features.m_transverse[i]);
        features.m_energy[i] = features.NormaliseEnergy(features.m_energy[i]);
        features.m_hitWidth[i] = features.NormaliseHitWidth(features.m_hitWidth[i]);
        features.m_theta[i] = features.NormaliseTheta(features.m_theta[i]);
        features.m_secVertex[i] = features.NormaliseSecVertex(features.m_secVertex[i]);
        features.m_eventHitSep[i] = features.NormaliseEventHitSep(features.m_eventHitSep[i]);
        features.m_clusterHitSep[i] = features.NormaliseClusterHitSep(features.m_clusterHitSep[i]);
        features.m_gapSep[i] = features.NormaliseGapSep(features.m_gapSep[i]);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DLClusterSplittingAlgorithm::GetAngle(const CartesianVector &position2D, const TwoDSlidingFitResult &clusterFit) const
{
    float thisL(0.f), thisT(0.f);
    clusterFit.GetLocalPosition(position2D, thisL, thisT);

    CartesianVector thisDirection(0.f, 0.f, 0.f);
    if (clusterFit.GetGlobalFitDirection(thisL, thisDirection) != STATUS_CODE_SUCCESS)
        return -4.f;

    const float openingAngleZ = CartesianVector(0.f, 0.f, 1.f).GetOpeningAngle(thisDirection);
    float openingAngleX = CartesianVector(1.f, 0.f, 0.f).GetOpeningAngle(thisDirection);
    openingAngleX *= (openingAngleZ > (M_PI * 0.5f)) ? (-1.f) : 1.f;

    return openingAngleX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DLClusterSplittingAlgorithm::GetDistanceToGap(const CartesianVector &position2D, const HitType hitType) const
{
    const DetectorGapList detectorGapList(this->GetPandora().GetGeometry()->GetDetectorGapList());

    if (detectorGapList.empty())
        return -1.f;

    float minDist(std::numeric_limits<float>::max());
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

        if (((hitType == TPC_VIEW_U) && (lineGapType == TPC_WIRE_GAP_VIEW_U)) ||
            ((hitType == TPC_VIEW_V) && (lineGapType == TPC_WIRE_GAP_VIEW_V)) ||
            ((hitType == TPC_VIEW_W) && (lineGapType == TPC_WIRE_GAP_VIEW_W)))
        {
            minDist = std::min(std::fabs(pLineGap->GetLineStartZ() - position2D.GetZ()), minDist);
            minDist = std::min(std::fabs(pLineGap->GetLineEndZ() - position2D.GetZ()), minDist);
        }
    }

    return minDist;
}

//------------------------------------------------------------------------------------------------------------------------------------------

    float DLClusterSplittingAlgorithm::GetDistanceToEventHit(HitKDTree2D &kdTree, const CaloHit *const pCaloHit, const CaloHitList &clusterPathHits) const
{
     // Collect close hits
    HitKDNode2DList foundHits;
    KDTreeBox searchRegionHits(build_2d_kd_search_region(pCaloHit, m_searchRegion1D, m_searchRegion1D));
    kdTree.search(searchRegionHits, foundHits);
    // Filter
    bool found(false);
    float minDistSq(std::numeric_limits<float>::max());
    for (const auto &hit : foundHits)
    {
        const CaloHit *const pFoundHit(hit.data);

        if (std::find(clusterPathHits.begin(), clusterPathHits.end(), pFoundHit) != clusterPathHits.end())
            continue;

        found = true;
        minDistSq = std::min(minDistSq, (pCaloHit->GetPositionVector() - pFoundHit->GetPositionVector()).GetMagnitudeSquared());
    }

    return (found ? std::sqrt(minDistSq) : -1.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DLClusterSplittingAlgorithm::GetDistanceToClusterHit(const ClusterPath::const_iterator &currentHit, 
    const ClusterPath::const_iterator &endIter) const
{
    const ClusterPath::const_iterator nextHit(std::next(currentHit));

    if (nextHit == endIter)
        return -1.f;

    return (currentHit->second.first->GetPositionVector() - nextHit->second.first->GetPositionVector()).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DLClusterSplittingAlgorithm::GetDistanceToSecVertex(const CaloHit *const pCaloHit, const HitType hitType) const
{
    if (!m_pSecVertexList)
        return -1.f;

    float bestSepSq(std::numeric_limits<float>::max());
    for (const Vertex *const pSecVertex : *m_pSecVertexList)
    {
        const CartesianVector secVtxPos(LArGeometryHelper::ProjectPosition(this->GetPandora(), pSecVertex->GetPosition(), hitType));
        bestSepSq = std::min(bestSepSq, (secVtxPos - pCaloHit->GetPositionVector()).GetMagnitudeSquared());
    }

    return std::sqrt(bestSepSq);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLClusterSplittingAlgorithm::GetWindows(Features &features, IntVector &splitIndices)
{
    const int m_windowLength(48);
    int sequenceLength(features.m_transverse.size());
    FloatVector windowStart, windowEnd;

    // If too small, then pad
    if (sequenceLength < m_windowLength)
    {
        for (int i = 0; i < (m_windowLength - sequenceLength); ++i)
        {
            features.m_transverse.push_back(-9999.9f);
            features.m_energy.push_back(-9999.9f);
            features.m_hitWidth.push_back(-9999.9f);
            features.m_theta.push_back(-9999.9f);
            features.m_secVertex.push_back(-9999.9f);
            features.m_gapSep.push_back(-9999.9f);
            features.m_eventHitSep.push_back(-9999.9f);
            features.m_clusterHitSep.push_back(-9999.9f);
        }

        windowStart.push_back(0);
        windowEnd.push_back(m_windowLength);
    }
    // If too big then split
    else if (sequenceLength > m_windowLength)
    {
        const int nWindows(std::floor(sequenceLength) / m_windowLength);

        for (int i = 0; i < nWindows; ++i)
        {
            windowStart.push_back(m_windowLength * i);
            windowEnd.push_back(m_windowLength * (i + 1));
        }

        if (sequenceLength % m_windowLength != 0)
        {
            windowStart.push_back(sequenceLength - m_windowLength);
            windowEnd.push_back(sequenceLength);
        }
    }

    // std::cout << "-------------------------------" << std::endl;
    // std::cout << "PRINTING WINDOWS" << std::endl;
    // std::cout << "sequence length: " << sequenceLength << std::endl;
    // std::cout << "-------------------------------" << std::endl;

    // (n_windows, sequence_length, n_features)
    LArDLHelper::TorchInput input;
    LArDLHelper::InitialiseInput({static_cast<int>(windowStart.size()), m_windowLength, 8}, input); 

    for (unsigned int i = 0; i < windowStart.size(); ++i)
    {
        for (int j = 0; j < m_windowLength; ++j)
        {
            input[i][j][0] = features.m_transverse.at(windowStart.at(i) + j);
            input[i][j][1] = features.m_energy.at(windowStart.at(i) + j);
            input[i][j][2] = features.m_hitWidth.at(windowStart.at(i) + j);
            input[i][j][3] = features.m_theta.at(windowStart.at(i) + j);
            input[i][j][4] = features.m_secVertex.at(windowStart.at(i) + j);
            input[i][j][5] = features.m_gapSep.at(windowStart.at(i) + j);
            input[i][j][6] = features.m_eventHitSep.at(windowStart.at(i) + j);
            input[i][j][7] = features.m_clusterHitSep.at(windowStart.at(i) + j);
        }
    }

    LArDLHelper::TorchOutput windowOutput, splitPosOutput;
    LArDLHelper::Forward(m_windowModel, {input}, windowOutput);
    LArDLHelper::Forward(m_splitPosModel, {input}, splitPosOutput);
    torch::TensorAccessor<float, 2> windowOutputAccessor = windowOutput.accessor<float, 2>();
    torch::TensorAccessor<float, 3> splitPosOutputAccessor = splitPosOutput.accessor<float, 3>();

    IntVector splitIndices_temp;
    FloatVector splitScores_temp;

    for (unsigned int i = 0; i < windowStart.size(); ++i)
    {
        // std::cout << "start: " << windowStart.at(i) << std::endl;
        // std::cout << "end: " << windowEnd.at(i) << std::endl;
        // std::cout << "not contaminated: " << static_cast<float>(windowOutputAccessor[i][0]) << std::endl;
        // std::cout << "contaminated: " << static_cast<float>(windowOutputAccessor[i][1]) << std::endl;
        // std::cout << "shower: " << static_cast<float>(windowOutputAccessor[i][2]) << std::endl;

        // Is contaminated?
        if (windowOutputAccessor[i][1] > 0.5)
        {
            // Is there a split point?
            for (int j = 0; j < m_windowLength; ++j)
            {
                const int sequenceIndex(windowStart.at(i) + j);
                const int nextWindowStart((i+1) == windowStart.size() ? std::numeric_limits<int>::max() : windowStart.at(i+1));

                // make sure we only use the last window for the overlap region...
                if (sequenceIndex >= nextWindowStart)
                    continue;

                if (static_cast<float>(splitPosOutputAccessor[i][j][0] > 0.5))
                {
                    //splitIndices.push_back(sequenceIndex);
                    splitIndices_temp.push_back(sequenceIndex);
                    splitScores_temp.push_back(static_cast<float>(splitPosOutputAccessor[i][j][0]));
                }
            }
        }
    }

    // The model will often identify consecutive splitting positions around the truth, so identify one point
    int previousIndex(-2), bestIndex(-1);
    float bestScore(-1.f);

    for (unsigned int iSplit = 0; iSplit < splitIndices_temp.size(); ++iSplit)
    {
        int thisIndex(splitIndices_temp.at(iSplit));
        float thisScore(splitScores_temp.at(iSplit));

        // If this is the first entry
        if (thisIndex == splitIndices_temp.front())
        {
            bestIndex = thisIndex;
            bestScore = thisScore;
        }
        else if (thisIndex == (previousIndex + 1))
        {
            if (thisScore > bestScore)
            {
                bestScore = thisScore;
                bestIndex = thisIndex;
            }
        }
        else
        {
            splitIndices.push_back(bestIndex);
            bestScore = -1.f;
            bestIndex = -1;
        }

        previousIndex = thisIndex;

        // If this is the last entry
        if (previousIndex == splitIndices_temp.back())
            splitIndices.push_back(bestIndex);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLClusterSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, 
        XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, 
        XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "SecVertexListName", m_secVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MinClusterHits", m_minClusterHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "SlidingWindow", m_slidingWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "LBinSize", m_lBinSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "SearchRegion1D", m_searchRegion1D));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "WindowModelName", m_windowModelName));
    m_windowModelName = LArFileHelper::FindFileInPath(m_windowModelName, "FW_SEARCH_PATH");
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_windowModelName, m_windowModel));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SplitPosModelName", m_splitPosModelName));
    m_splitPosModelName = LArFileHelper::FindFileInPath(m_splitPosModelName, "FW_SEARCH_PATH");
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_splitPosModelName, m_splitPosModel));


    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
