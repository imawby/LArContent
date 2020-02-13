/**
 *  @file   HitWidthClusterMergingAlgorithm.h
 *
 *  @brief  Header file for the hit width cluster merging algorithm class.
 *
 *  $Log: $
 */

#ifndef LAR_HIT_WIDTH_CLUSTER_MERGING_ALGORITHM_H
#define LAR_HIT_WIDTH_CLUSTER_MERGING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/ClusterAssociationAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

namespace lar_content
{


/**
 *  @brief HitWidthClusterMergingAlgorithm class
 */
  class HitWidthClusterMergingAlgorithm : public pandora::Algorithm//public ClusterAssociationAlgorithm
{

public:

    /**
     *  @brief  Default constructor
     */
    HitWidthClusterMergingAlgorithm();


    class ClusterFit
    {
    public:
      /**
      * @brief  Default constructor
      */
      ClusterFit(const pandora::Cluster *const pCluster);

      pandora::CartesianVector m_lowerXExtrema;
      pandora::CartesianVector m_higherXExtrema;
      
      struct ClusterPositionSort{
	  ClusterPositionSort(const pandora::CartesianVector referencePoint) : m_referencePoint(referencePoint) {}
          bool operator() (const pandora::CartesianVector &lhs, const pandora::CartesianVector &rhs) {
              return (m_referencePoint.GetDistanceSquared(lhs) < m_referencePoint.GetDistanceSquared(rhs));
	  } 
	  const pandora::CartesianVector m_referencePoint;
      };


      typedef std::multimap<pandora::CartesianVector, float, ClusterPositionSort> ClusterPositionToWeightMap; 
      

      unsigned int m_numCaloHits;
      float m_totalWeight;

      ClusterPositionSort m_currentClusterSort;
      ClusterPositionSort m_testClusterSort;


      ClusterPositionToWeightMap m_currentClusterPositionToWeightMap;
      ClusterPositionToWeightMap m_testClusterPositionToWeightMap;


    };





private:
  
  //void GetWeightedSubGradient(const pandora::Cluster *const pCluster, bool isTransverse, bool isCurrent, unsigned int numFittingPoints, pandora::CartesianVector &direction, float &chiSquared) const;



  pandora::StatusCode Run();
  
  //void TestGetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;

  //void TestPopulateClusterAssociationMap(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const;

  //bool TestAreClustersAssociated(const pandora::Cluster *const pCurrentCluster, const pandora::Cluster *const pTestCluster) const;
  
  void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;

  //void PopulateClusterAssociationMap(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const;

  //bool IsExtremalCluster(const bool isForward, const pandora::Cluster *const pCurrentCluster,  const pandora::Cluster *const pTestCluster) const;

  //bool AreClustersAssociated(const pandora::Cluster *const pCurrentCluster, const pandora::Cluster *const pTestCluster) const;

  //void GetWeightedGradient(const pandora::Cluster *const pCluster, bool isTransverse, pandora::CartesianVector &direction, pandora::CartesianVector &intercept, float &chiSquared) const;

  //pandora::CartesianVector GetClusterDirection(const pandora::Cluster *const pCluster) const;

  //pandora::CartesianVector GetClusterZIntercept(const pandora::Cluster *const pCluster) const;
 
  
  //THESE REALLY BELONG IN THE HELPER CLASS

  static pandora::CartesianVector GetExtremalCoordinatesLowerX(const pandora::Cluster *const pCluster);
  static pandora::CartesianVector GetExtremalCoordinatesHigherX(const pandora::Cluster *const pCluster) ;
  
  static void GetExtremalCoordinates(const pandora::Cluster *const pCluster, pandora::CartesianVector &lowerXCoordinate, pandora::CartesianVector &higherXCoordinate);

  static void GetExtremalCoordinates(const pandora::OrderedCaloHitList &orderedCaloHitList, pandora::CartesianVector &lowerXCoordinate, pandora::CartesianVector &higherXCoordinate);

  static void GetExtremalCoordinates(const pandora::CartesianPointVector &coordinateVector, pandora::CartesianVector &lowerXCoordinate, pandora::CartesianVector &higherXCoordinate);

  static void GetExtremalCoordinatesX(const pandora::Cluster *const pCluster, float &minX, float &maxX);

  static bool SortByX(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

  static float GetMinX(const pandora::Cluster *const pCluster);

  static bool SortByMaxX(const pandora::Cluster *const pLhs, const pandora::Cluster *const pRhs);

  static float GetTotalClusterWeight(const pandora::Cluster *const pCluster);
  
  //void GetPositionToWeightMap(const pandora::Cluster *const pCluster, LArClusterHelper::PositionToWeightMap &positionToWeightMap);

  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

 
  std::string m_clusterListName;
  float m_minClusterWeight;
  float m_maxXMergeDistance; //Distance either side of point
  float m_maxZMergeDistance; //Distance either side of point
  float m_maxMergeCosOpeningAngle; 

  //CartesianVector m_sortReferencePoint;

};

} //namespace lar_content

#endif //LAR_HIT_WIDTH_CLUSTER_MERGING_ALGORITHM_H
