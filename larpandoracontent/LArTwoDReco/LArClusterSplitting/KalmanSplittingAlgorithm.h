/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanSplittingAlgorithm.h
 *
 *  @brief  Header file for the two dimensional sliding fit splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_KALMAN_SPLITTING_ALGORITHM_H
#define LAR_KALMAN_SPLITTING_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/ClusterSplittingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  KalmanSplittingAlgorithm class
 */
class KalmanSplittingAlgorithm : public ClusterSplittingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    KalmanSplittingAlgorithm();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StatusCode DivideCaloHits(
        const pandora::Cluster *const pCluster, pandora::CaloHitList &firstCaloHitList, pandora::CaloHitList &secondCaloHitList) const;

    pandora::StatusCode FindBestSplitPosition(const pandora::Cluster *const pCluster, pandora::CartesianVector &splitPosition) const;

    /* pandora::StatusCode DivideCaloHits(const TwoDSlidingFitResult &slidingFitResult, const pandora::CartesianVector &splitPosition, */
    /*     pandora::CaloHitList &firstCaloHitList, pandora::CaloHitList &secondCaloHitList) const; */

    float m_minClusterLength;
};

} // namespace lar_content

#endif // #ifndef LAR_KALMAN_SPLITTING_ALGORITHM_H
