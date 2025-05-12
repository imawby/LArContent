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

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanFit.h"
#include "larpandoracontent/LArUtility/KalmanFilter.h"

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

    pandora::StatusCode DivideCaloHits(const pandora::Cluster *const pCluster, pandora::CaloHitList &firstCaloHitList, pandora::CaloHitList &secondCaloHitList) const;

    bool IsTargetCluster(const pandora::Cluster *const pCluster) const;

    pandora::StatusCode GetFitSeed(const std::map<int, pandora::CaloHitList> &caloHitWireMap, Eigen::VectorXd &init, int &startWireID) const;

    KalmanFit PerformKalmanFit(const pandora::Cluster *const pCluster) const;

    pandora::StatusCode FindMatchedClusterPosition(const pandora::CaloHitList &caloHitList, const KalmanFit &kalmanFit,
        pandora::CartesianPointVector &matchedPositions) const ;

    pandora::StatusCode FindBestSplitPosition(const pandora::Cluster *const pCluster, const KalmanFit &kalmanFit, 
        pandora::CartesianVector &splitPosition) const;

    void FollowRoute(const std::map<int, pandora::CaloHitList> &caloHitWireMap, KalmanFit &kalmanFit) const;

    pandora::CartesianVector GetMedian(pandora::CartesianPointVector &cartesianPointVector) const;

    float GetSTD(const pandora::CartesianPointVector &kalmanFit_dir) const;

    pandora::StatusCode DivideCaloHits(const TwoDSlidingFitResult &slidingFitResult, const pandora::CartesianVector &splitPosition, 
        pandora::CaloHitList &firstCaloHitList, pandora::CaloHitList &secondCaloHitList) const;

    float m_minClusterLength;
    float m_kalmanDelta;
    float m_kalmanProcessVarCoeff;
    float m_kalmanMeasurementVarCoeff;
    float m_minTransSeparation;
    int m_segmentWindows;
    float m_minDeviation;
    float m_maxSpread;
};

} // namespace lar_content

#endif // #ifndef LAR_KALMAN_SPLITTING_ALGORITHM_H
