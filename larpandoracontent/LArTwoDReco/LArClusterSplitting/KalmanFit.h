/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanFilter.h
 *
 *  @brief  Kalman fit
 *
 *  $Log: $
 */
#ifndef LAR_KALMAN_FIT_H
#define LAR_KALMAN_FIT_H 1

#include "Api/PandoraApi.h"
#include "larpandoracontent/LArUtility/KalmanFilter.h"

namespace lar_content
{
  class KalmanFit
  {

  public:
      KalmanFit(const KalmanFilter2D &kalmanFilter, const int currentWireID);
      void SaveStep(const int currentWireID);
      void AddPositionAndUpdate(const pandora::CartesianVector &position);

      KalmanFilter2D m_kalmanFilter;
      pandora::CartesianPointVector m_positions;
      pandora::CartesianPointVector m_directions;
      int m_currentWireID; // last wireID assessed
  };

  typedef std::vector<KalmanFit> KalmanFitVector;
}

#endif // #ifndef LAR_KALMAN_FIT_H
