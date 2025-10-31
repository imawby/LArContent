/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the two dimensional sliding fit splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArUtility/KalmanFilter.h"
#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanSplittingAlgorithm.h"



using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

KalmanSplittingAlgorithm::KalmanSplittingAlgorithm() :
    m_secVertexListName("SecondaryVertices3D"),
    m_pSecVertexList(nullptr),
    m_minClusterHits(50),
    m_slidingWindow(20),
    m_lBinSize(0.5f),
    m_searchRegion1D(20.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanSplittingAlgorithm::Run()
{
    // Get view hits
    const CaloHitList *pCaloHitList(nullptr);
    if (PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_SUCCESS;

    if ((!pCaloHitList) || pCaloHitList->empty())
        return STATUS_CODE_SUCCESS;

    // Get secondary vertices (it's okay if the list is empty)
    PandoraContentApi::GetList(*this, m_secVertexListName, m_pSecVertexList);

    // Now we've setup, run base alg
    StatusCode stat(ClusterSplittingAlgorithm::Run());

    if (stat != STATUS_CODE_SUCCESS)
    {
        m_pSecVertexList = nullptr;
        return stat;
    }

    m_pSecVertexList = nullptr;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanSplittingAlgorithm::DivideCaloHits(const Cluster *const pCluster, CaloHitList &/*firstHitList*/, 
    CaloHitList &/*secondHitList*/) const
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

        // Get windows
        this->GetWindows(features);
    }
    catch (...)
    {
        return STATUS_CODE_NOT_FOUND;
    }

    // Do not want to run twice?
    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanSplittingAlgorithm::FindPath(const Cluster *const pCluster, const TwoDSlidingFitResult &clusterFit, ClusterPath &clusterPath) const
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

void KalmanSplittingAlgorithm::FillFeatures(const ClusterPath &clusterPath, const TwoDSlidingFitResult &clusterFit, Features &features) const
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

    // Fill features
    for (ClusterPath::const_iterator iter = clusterPath.begin(); iter != clusterPath.end(); ++iter)
    {
        if (!processedFirst)
        {
            features.m_theta.push_back(features.NormaliseTheta(-4.f));
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
            features.m_theta.push_back(features.NormaliseTheta(openingAngleL));
        }

        const CaloHit *const pPathCaloHit(iter->second.first);
        const float energyFrac(pPathCaloHit->GetElectromagneticEnergy() / totalEnergy);

        features.m_transverse.push_back(
            features.NormaliseTransverse(iter->second.second));
        features.m_energy.push_back(
            features.NormaliseEnergy(features.m_energy.empty() ? energyFrac : (energyFrac + features.m_energy.back())));
        features.m_hitWidth.push_back(
            features.NormaliseHitWidth(pPathCaloHit->GetCellSize1()));
        features.m_secVertex.push_back(
            features.NormaliseSecVertex(this->GetDistanceToSecVertex(pPathCaloHit, pPathCaloHit->GetHitType())));
        features.m_eventHitSep.push_back(
            features.NormaliseEventHitSep(this->GetDistanceToEventHit(kdTree, pPathCaloHit, clusterPathHits)));
        features.m_clusterHitSep.push_back(
            features.NormaliseClusterHitSep(this->GetDistanceToClusterHit(iter, clusterPath.end())));
        features.m_gapSep.push_back(
            features.NormaliseGapSep(this->GetDistanceToGap(pPathCaloHit->GetPositionVector(), pPathCaloHit->GetHitType())));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float KalmanSplittingAlgorithm::GetAngle(const CartesianVector &position2D, const TwoDSlidingFitResult &clusterFit) const
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

float KalmanSplittingAlgorithm::GetDistanceToGap(const CartesianVector &position2D, const HitType hitType) const
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

    float KalmanSplittingAlgorithm::GetDistanceToEventHit(HitKDTree2D &kdTree, const CaloHit *const pCaloHit, const CaloHitList &clusterPathHits) const
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

float KalmanSplittingAlgorithm::GetDistanceToClusterHit(const ClusterPath::const_iterator &currentHit, 
    const ClusterPath::const_iterator &endIter) const
{
    const ClusterPath::const_iterator nextHit(std::next(currentHit));

    if (nextHit == endIter)
        return -1.f;

    return (currentHit->second.first->GetPositionVector() - nextHit->second.first->GetPositionVector()).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float KalmanSplittingAlgorithm::GetDistanceToSecVertex(const CaloHit *const pCaloHit, const HitType hitType) const
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

void KalmanSplittingAlgorithm::GetWindows(Features &features) const
{
    const int m_windowLength(48);
    int sequenceLength(features.m_transverse.size());
    FloatVector windowStart, windowEnd;
    // std::vector<FloatVector> positions;
    // FloatVector windowScores;
    // std::vector<FloatVector> splitScores;

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

    std::cout << "-------------------------------" << std::endl;
    std::cout << "PRINTING WINDOWS" << std::endl;
    std::cout << "sequence length: " << sequenceLength << std::endl;
    std::cout << "-------------------------------" << std::endl;

    for (unsigned int i = 0; i < windowStart.size(); ++i)
    {
        std::cout << "start: " << windowStart.at(i) << std::endl;
        std::cout << "end: " << windowEnd.at(i) << std::endl;

        FloatVector v2 = FloatVector(features.m_transverse.begin() + windowStart.at(i), features.m_transverse.begin() + windowEnd.at(i));
        std::cout << "length: " << v2.size() << std::endl;

    }


    // LArDLHelper::TorchInput input;
    // LArDLHelper::InitialiseInput({1, 6}, input);

    // int insertIndex(0);

    // for (const FloatVector &edgeOutput : {outputUp, outputDown})
    // {
    //     for (int i = 0; i < 3; ++i)
    //     {
    //         input[0][insertIndex] = edgeOutput.at(i);
    //         ++insertIndex;
    //     }
    // }

    // LArDLHelper::TorchOutput output;
    // LArDLHelper::Forward(m_primaryTrackClassifierModel, {input}, output);
    // torch::TensorAccessor<float, 2> outputAccessor = output.accessor<float, 2>();




}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, 
        XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

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

    return ClusterSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
