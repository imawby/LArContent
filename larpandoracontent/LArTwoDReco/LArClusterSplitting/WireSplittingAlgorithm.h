/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/WireSplittingAlgorithm.h
 *
 *  @brief  Header file for the wire splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_WIRE_SPLITTING_ALGORITHM_H
#define LAR_WIRE_SPLITTING_ALGORITHM_H 1

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/ClusterSplittingAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  WireSplittingAlgorithm class
 */
class WireSplittingAlgorithm : public ClusterSplittingAlgorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    WireSplittingAlgorithm();

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StatusCode DivideCaloHits(
        const pandora::Cluster *const pCluster, pandora::CaloHitList &firstCaloHitList, pandora::CaloHitList &secondCaloHitList) const;

    float m_minClusterLength;            ///<
};

} // namespace lar_content

#endif // #ifndef LAR_WIRE_SPLITTING_ALGORITHM_H
