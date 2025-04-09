/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the two dimensional sliding fit splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

KalmanSplittingAlgorithm::KalmanSplittingAlgorithm() :
    m_minClusterLength(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanSplittingAlgorithm::DivideCaloHits(const Cluster *const pCluster, CaloHitList &firstHitList, CaloHitList &secondHitList) const
{
    if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
        return STATUS_CODE_NOT_FOUND;


    CartesianVector splitPosition(0.f, 0.f, 0.f);
    this->FindBestSplitPosition(pCluster, splitPosition);


    // try
    // {
    //     if (STATUS_CODE_SUCCESS == this->FindBestSplitPosition(pCluster, splitPosition))
    //     {
    //         return this->DivideCaloHits(slidingFitResult, splitPosition, firstHitList, secondHitList);
    //     }
    // }
    // catch (StatusCodeException &statusCodeException)
    // {
    //     if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
    //         throw statusCodeException;
    // }

    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanSplittingAlgorithm::FindBestSplitPosition(const Cluster *const pCluster, CartesianVector &/*splitPosition*/) const
{
    // Get CaloHits
    CaloHitList caloHitList;
    LArClusterHelper::GetAllHits(pCluster, caloHitList);

    // Begin with a PCA
    CartesianVector centroid(0.f, 0.f, 0.f);
    LArPcaHelper::EigenVectors eigenVecs;
    LArPcaHelper::EigenValues eigenValues(0.f, 0.f, 0.f);
    LArPcaHelper::RunPca(caloHitList, centroid, eigenValues, eigenVecs);

    // Make sure sensible?

    // Order hits wrt l
    CaloHitVector caloHitVector;
    for (const CaloHit *const pCaloHit : caloHitList)
        caloHitVector.emplace_back(pCaloHit);

    std::sort(caloHitVector.begin(), caloHitVector.end(), [&](const CaloHit *const i, const CaloHit *const j) 
        { 
            const float l_i(eigenVecs.front().GetDotProduct(i->GetPositionVector() - centroid));
            const float l_j(eigenVecs.front().GetDotProduct(j->GetPositionVector() - centroid));

            return l_i < l_j;
        });

    // Now build kalman filter
    const float m_kalmanDelta(1.f);
    const float m_kalmanProcessVarCoeff(1.f);
    const float m_kalmanMeasurementVarCoeff(1.f);
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const HitType view(caloHitVector.front()->GetHitType());
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
    const float processVariance{m_kalmanProcessVarCoeff * pitch * pitch};
    const float measurementVariance{m_kalmanMeasurementVarCoeff * pitch * pitch};

    KalmanFilter kalman(m_kalmanDelta, processVariance, measurementVariance, caloHitVector.front());

    kalman.Predict();


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

// StatusCode KalmanSplittingAlgorithm::DivideCaloHits(const TwoDSlidingFitResult &slidingFitResult,
//     const CartesianVector &splitPosition, CaloHitList &firstCaloHitList, CaloHitList &secondCaloHitList) const
// {
//     return STATUS_CODE_SUCCESS;
// }

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", m_minClusterLength));

    return ClusterSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
