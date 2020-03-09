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
  class HitWidthClusterMergingAlgorithm : public ClusterAssociationAlgorithm
{

public:

    /**
     *  @brief  Default constructor
     */
    HitWidthClusterMergingAlgorithm();

    
    class ClusterFitParameters
    {
    public:
      /**
      * @brief  Default constructor
      */
      ClusterFitParameters(const pandora::Cluster *const pCluster, const float maxConsituentHitWidth);
      
      struct ClusterPositionSort{

          ClusterPositionSort(const pandora::CartesianVector referencePoint) : m_referencePoint(referencePoint) {}

          bool operator() (const pandora::CartesianVector &lhs, const pandora::CartesianVector &rhs) {
              return (m_referencePoint.GetDistanceSquared(lhs) < m_referencePoint.GetDistanceSquared(rhs));
	  } 

	  const pandora::CartesianVector m_referencePoint;
      };


      typedef std::multimap<pandora::CartesianVector, float, ClusterPositionSort> ClusterPositionToWeightMap; 
      
      const pandora::Cluster *m_pCluster;                 ///< The address of the cluster

      unsigned int m_numCaloHits;
      float m_totalWeight;

      pandora::CartesianVector m_lowerXExtrema;
      pandora::CartesianVector m_higherXExtrema;

      ClusterPositionSort m_currentClusterSort;
      ClusterPositionSort m_testClusterSort;

      ClusterPositionToWeightMap m_currentClusterPositionToWeightMap;
      ClusterPositionToWeightMap m_testClusterPositionToWeightMap;



    };

    typedef std::unordered_map<const pandora::Cluster*, const ClusterFitParameters> ClusterToFitParametersMap;

private:
  
    //pandora::StatusCode Run();
  
  
  void GetListOfCleanClusters(const pandora::ClusterList *const pClusterList, pandora::ClusterVector &clusterVector) const;

  void PopulateClusterAssociationMap(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const;

  void CleanupClusterAssociations(const pandora::ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const;

  bool AreClustersAssociated(const ClusterFitParameters &currentFitParameters, const ClusterFitParameters &testFitParameters) const;

  bool IsExtremalCluster(const bool isForward, const pandora::Cluster *const pCurrentCluster,  const pandora::Cluster *const pTestCluster) const;

  pandora::CartesianVector GetClusterDirection(const ClusterFitParameters::ClusterPositionToWeightMap &clusterPositionToWeightMap, unsigned int clusterCaloHits) const;

  pandora::CartesianVector GetClusterZIntercept(const ClusterFitParameters::ClusterPositionToWeightMap &clusterPositionToWeightMap, unsigned int clusterCaloHits) const;

  void GetWeightedGradient(const ClusterFitParameters::ClusterPositionToWeightMap &clusterPositionToWeightMap, bool isTransverse, pandora::CartesianVector &direction, pandora::CartesianVector &intercept, float &chiSquared, float clusterCaloHits) const;


  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


  std::string m_clusterListName;
  float m_minClusterWeight;
  float m_maxXMergeDistance; //Distance either side of point
  float m_maxZMergeDistance; //Distance either side of point
  float m_maxMergeCosOpeningAngle; 

  float m_maxConstituentHitWidth;

  bool m_fitToFullCluster;
  float m_fittingSampleWeight;

  ClusterToFitParametersMap m_clusterToFitParametersMap;
};

} //namespace lar_content

#endif //LAR_HIT_WIDTH_CLUSTER_MERGING_ALGORITHM_H



