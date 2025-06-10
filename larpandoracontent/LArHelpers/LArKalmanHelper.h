/**
 *  @file   larpandoracontent/LArHelpers/LArKalmanHelper.h
 *
 *  @brief  Header file for the kalman helper class.
 *
 *  $Log: $
 */
#ifndef LAR_KALMAN_HELPER_H
#define LAR_KALMAN_HELPER_H 1

#include "Eigen/Dense"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanFit.h"

#include "Objects/Cluster.h"

namespace lar_content
{

/**
 *  @brief  LArKalmanHelper class
 */
class LArKalmanHelper
{
public:
    static pandora::StatusCode GetFitSeed(const std::map<int, pandora::CaloHitList> &caloHitWireMap, Eigen::VectorXd &init, int &startWireID);

    static pandora::StatusCode FindMatchedClusterPosition(const pandora::CaloHitList &caloHitList, const KalmanFit &kalmanFit,
        const float minTransSep, pandora::CartesianPointVector &matchedPositions);

    static void FollowRoute(const std::map<int, pandora::CaloHitList> &caloHitWireMap, const float minTransSep, KalmanFit &kalmanFit);

    static float GetSTD(const pandora::CartesianPointVector &directionVector);

    static KalmanFit PerformKalmanFit(const pandora::Pandora &pandora, const pandora::Cluster *const pCluster, const float minTransSep);

};

} // namespace lar_content

#endif // #ifndef LAR_KALMAN_HELPER_H
