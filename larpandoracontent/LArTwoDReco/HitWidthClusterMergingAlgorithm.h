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

namespace lar_content
{

/**
 *  @brief HitWidthClusterMergingAlgorithm class
 */
class HitWidthClusterMergingAlgorithm : public pandora::Algorithm
{

public:
    /**
     *  @brief  Default constructor
     */
  HitWidthClusterMergingAlgorithm();

private:

  pandora::StatusCode Run();

  void GetExtremalCoordinates(const pandora::Cluster *const pCluster, pandora::CartesianVector &innerCoordinate, pandora::CartesianVector &outerCoordinate);

  void GetExtremalCoordinates(const pandora::OrderedCaloHitList &orderedCaloHitList, pandora::CartesianVector &innerCoordinate, pandora::CartesianVector &outerCoordinate);

  pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

};

} //namespace lar_content

#endif //LAR_HIT_WIDTH_CLUSTER_MERGING_ALGORITHM_H
