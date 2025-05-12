/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanFit.cc
 *
 *  @brief  Kalman fit object
 *
 *  $Log: $
 */

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanFit.h"

using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

KalmanFit::KalmanFit(const KalmanFilter2D &kalmanFilter, const int currentWireID) :
    m_kalmanFilter(kalmanFilter),
    m_positions(CartesianPointVector()),
    m_directions(CartesianPointVector()),
    m_currentWireID(currentWireID)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanFit::SaveStep(const int currentWireID)
{
    const CartesianVector kalmanPos(m_kalmanFilter.GetPosition()(0), 0.f, m_kalmanFilter.GetPosition()(1));
    const CartesianVector kalmanDir(m_kalmanFilter.GetDirection()(0), 0.f, m_kalmanFilter.GetDirection()(1));

    m_positions.push_back(kalmanPos);
    m_directions.push_back(kalmanDir);
    m_currentWireID = currentWireID;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanFit::AddPositionAndUpdate(const CartesianVector &position)
{
    Eigen::VectorXd eigenXd(2);
    eigenXd << position.GetX(), position.GetZ();
    m_kalmanFilter.Update(eigenXd);
}

//------------------------------------------------------------------------------------------------------------------------------------------

}
