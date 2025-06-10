/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the two dimensional sliding fit splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"
#include "larpandoracontent/LArHelpers/LArKalmanHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArUtility/KalmanFilter.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanFit.h"
#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

KalmanSplittingAlgorithm::KalmanSplittingAlgorithm() :
    m_minClusterLength(10.f),
    m_kalmanDelta(1.f),
    m_kalmanProcessVarCoeff(1.f),
    m_kalmanMeasurementVarCoeff(1.f),
    m_minTransSeparation(5.f),
    m_segmentWindows(6),
    m_minDeviation(10.f),
    m_maxSpread(0.1f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanSplittingAlgorithm::DivideCaloHits(const Cluster *const pCluster, CaloHitList &firstHitList, CaloHitList &secondHitList) const
{
    // Do we want to work on this cluster?
    if (!this->IsTargetCluster(pCluster))
        return STATUS_CODE_NOT_FOUND;

    //////////////////////////////////
    // PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    //////////////////////////////////

    // Get KalmanFit
    try
    {
        //////////////////////////////////
        // ClusterList visualiseClusters({pCluster});
        // PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &visualiseClusters, "Cluster", (pCluster->GetParticleId() == MU_MINUS ? BLUE : RED)));
        // PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        //////////////////////////////////

        KalmanFit kalmanFit(LArKalmanHelper::PerformKalmanFit(this->GetPandora(), pCluster, m_minTransSeparation));

        if (kalmanFit.m_positions.empty())
            return STATUS_CODE_NOT_FOUND;

        // Do sliding fit
        const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
        const HitType view(LArClusterHelper::GetClusterHitType(pCluster));
        const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
        const TwoDSlidingFitResult twoDSlidingFitResult(pCluster, 20, pitch);

        // Search for split position
        CartesianVector splitPosition(0.f, 0.f, 0.f);

        if (STATUS_CODE_SUCCESS == this->FindBestSplitPosition(pCluster, kalmanFit, splitPosition))
        {
            return this->DivideCaloHits(twoDSlidingFitResult, splitPosition, firstHitList, secondHitList);
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
            throw statusCodeException;
    }

    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool KalmanSplittingAlgorithm::IsTargetCluster(const Cluster *const pCluster) const
{
    if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
        return false;

    const float m_maxConstituentHitWidth(0.5f), m_hitWidthScalingFactor(1.f), m_minClusterSparseness(0.5f);
    const unsigned int numberOfProposedConstituentHits(LArHitWidthHelper::GetNProposedConstituentHits(
        pCluster, m_maxConstituentHitWidth, m_hitWidthScalingFactor));

    if (numberOfProposedConstituentHits == 0)
        return false;

    // Avoid wide hit width clusters
    // clusterSparseness [0 -> 1] where a higher value indicates sparseness
    const float clusterSparseness(1.f - (static_cast<float>(pCluster->GetNCaloHits()) / static_cast<float>(numberOfProposedConstituentHits)));

    if (clusterSparseness > m_minClusterSparseness)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanSplittingAlgorithm::FindBestSplitPosition(const Cluster *const pCluster, const KalmanFit &kalmanFit, 
    CartesianVector &splitPosition) const
{
    // Probe kalman fit
    const int nFitPoints(kalmanFit.m_positions.size());
    if (nFitPoints < ((2 * m_segmentWindows) + 1)) { return STATUS_CODE_NOT_FOUND; };

    bool splitFound(false);
    float maxDeviation(std::numeric_limits<float>::min());

    for (unsigned int i = (m_segmentWindows + 1); i < static_cast<unsigned int>(nFitPoints - m_segmentWindows); ++i)
    {
        // Get before/after segments
        CartesianPointVector beforeSeg, afterSeg;
        for (int j = 1; j <= m_segmentWindows; j++)
        {
            beforeSeg.push_back(kalmanFit.m_directions.at(i - j)); 
            afterSeg.push_back(kalmanFit.m_directions.at(i + j));
        }

        // Get median direction from before/after region
        const CartesianVector beforeMedian(this->GetMedian(beforeSeg));
        const CartesianVector afterMedian(this->GetMedian(afterSeg));

        // Get splitting metrucs
        const float direct((kalmanFit.m_directions.at(i-1).GetOpeningAngle(kalmanFit.m_directions.at(i+1))) / 3.14 * 180.f);
        const float eitherSide(beforeMedian.GetOpeningAngle(afterMedian) / 3.14 * 180.f);
        const float afterSigma(LArKalmanHelper::GetSTD(afterSeg)); // make sure that after segment looks track-like...
        const float beforeSigma(LArKalmanHelper::GetSTD(beforeSeg)); // make sure that before segment looks track-like...

        // Measure angular deviation
        if ((eitherSide > m_minDeviation) && (beforeSigma < m_maxSpread) && (afterSigma < m_maxSpread) && (direct > m_minDeviation))
        { 
            if (eitherSide > maxDeviation)
            {
                splitFound = true;
                maxDeviation = eitherSide;
                splitPosition = kalmanFit.m_positions.at(i);
            }
        }
   } 

    // if (splitFound)
    // {
    //     PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &splitPosition, "splitPosition", BLUE, 2));
    //     PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    // }

    return splitFound ? STATUS_CODE_SUCCESS : STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector KalmanSplittingAlgorithm::GetMedian(CartesianPointVector &cartesianPointVector) const
{
    // Sort wrt x-axis theta
    std::sort(cartesianPointVector.begin(), cartesianPointVector.end(), [](const CartesianVector &lhs, const CartesianVector &rhs) 
        { 
            CartesianVector xAxis(1.f, 0.f, 0.f);
            float lhsAngle(xAxis.GetOpeningAngle(lhs));
            float rhsAngle(xAxis.GetOpeningAngle(rhs));

            if (lhs.GetY() < 0.f)
                lhsAngle += M_PI;

            if (rhs.GetY() < 0.f)
                rhsAngle += M_PI;

            return lhsAngle > rhsAngle;
        }
    );

    const int medianIndex = std::floor(cartesianPointVector.size() / 2);

    return cartesianPointVector.at(medianIndex);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanSplittingAlgorithm::DivideCaloHits(const TwoDSlidingFitResult &slidingFitResult,
    const CartesianVector &splitPosition, CaloHitList &firstCaloHitList, CaloHitList &secondCaloHitList) const
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

StatusCode KalmanSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "KalmanDelta", m_kalmanDelta));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "KalmanProcessVarCoeff", m_kalmanProcessVarCoeff));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "KalmanMeasurementVarCoeff", m_kalmanMeasurementVarCoeff));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinTransSeparation", m_minTransSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SegmentWindows", m_segmentWindows));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinDeviation", m_minDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSpread", m_maxSpread));

    return ClusterSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
