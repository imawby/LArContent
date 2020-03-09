/**
 *  @file   larpandoracontent/LArHelpers/LArHitWidthHelper.h
 *
 *  @brief  Header file for the lar hit width helper class.
 *
 *  $Log: $
 */

#ifndef LAR_HIT_WIDTH_HELPER_H
#define LAR_HIT_WIDTH_HELPER_H 1


#include "Objects/Cluster.h"

namespace lar_content
{
  

/**
 *  @brief  LArHitWidthHelper class
 */
class LArHitWidthHelper
{
public:

    static void GetExtremalCoordinatesX(const pandora::Cluster *const pCluster, pandora::CartesianVector &lowerXCoordinate, pandora::CartesianVector &higherXCoordinate, const float maxConstituentHitWidth);

    static pandora::CartesianPointVector GetConstituentHits(const pandora::Cluster *const pCluster, const float maxConstituentHitWidth);

    static pandora::CartesianPointVector GetUniformConstituentHits(const pandora::Cluster *const pCluster, const float constituentHitWidth);
 
    static pandora::CartesianVector GetExtremalCoordinatesLowerX(const pandora::Cluster *const pCluster, const float maxConstituentHitWidth);

    static pandora::CartesianVector GetExtremalCoordinatesHigherX(const pandora::Cluster *const pCluster, const float maxConstituentHitWidth) ;

  
  static bool SortByMaxX(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs,  const float maxConstituentHitWidth);


  static float GetTotalClusterWeight(const pandora::Cluster *const pCluster);
 


  struct MaxXPositionSort{

      MaxXPositionSort(const float maxConsituentHitWidth) : m_maxConstituentHitWidth(maxConsituentHitWidth) {}

      bool operator() (const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs) {
	  pandora::CartesianVector lhsLowerXEdge(0,0,0);
	  pandora::CartesianVector lhsUpperXEdge(0,0,0);

	  pandora::CartesianVector rhsLowerXEdge(0,0,0);
	  pandora::CartesianVector rhsUpperXEdge(0,0,0);

	  LArHitWidthHelper::GetExtremalCoordinatesX(pLhs, lhsLowerXEdge, lhsUpperXEdge, m_maxConstituentHitWidth);
	  LArHitWidthHelper::GetExtremalCoordinatesX(pRhs, rhsLowerXEdge, rhsUpperXEdge, m_maxConstituentHitWidth);
  
          return (lhsUpperXEdge.GetX() < rhsUpperXEdge.GetX());
      }

      const float m_maxConstituentHitWidth;
  };




};

} // namespace lar_content


#endif // #ifndef LAR_MC_PARTICLE_HELPER_H
