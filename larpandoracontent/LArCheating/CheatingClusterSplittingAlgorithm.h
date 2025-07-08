/**
 *  @file   larpandoracontent/LArCheating/CheatingClusterSplittingAlgorithm.h
 *
 *  @brief  Header file for the cheating cluster splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_CLUSTER_SPLITTING_ALGORITHM_H
#define LAR_CHEATING_CLUSTER_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  CheatingClusterSplittingAlgorithm class
 */
class CheatingClusterSplittingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingClusterSplittingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::unordered_map<const pandora::MCParticle *, pandora::CaloHitList> MCParticleToHitListMap;
    typedef std::unordered_map<const pandora::CaloHit *, const pandora::MCParticle *> HitToMCParticleMap;
    typedef std::unordered_map<const pandora::Cluster *, const pandora::MCParticle *> ClusterToMCParticleMap;
    typedef std::unordered_map<const pandora::Cluster *, pandora::MCParticleList> ClusterToMCParticleListMap;
    typedef std::unordered_map<const pandora::MCParticle *, pandora::CartesianVector> MCParticleSecVertexMap;

    void FillPandoraMaps(const pandora::ClusterList *const pClusterList, const pandora::CaloHitList *const pCaloHitList, 
                         const pandora::MCParticleList *const pMCParticleList, const pandora::VertexList *const pSecVertexList,
                         MCParticleToHitListMap &mcParticleToHitListMap, HitToMCParticleMap &hitToMCParticleMap, ClusterToMCParticleMap &clusterToMCParticleMap,
                         ClusterToMCParticleListMap &clusterToMCParticleListMap, MCParticleSecVertexMap &mcParticleSecVertexMap);

    bool SplitClusters(const pandora::ClusterList *const pClusterList, MCParticleSecVertexMap &mcParticleSecVertexMap, 
                       ClusterToMCParticleMap &clusterToMCParticleMap, ClusterToMCParticleListMap &clusterToMCParticleListMap,
                       MCParticleToHitListMap &mcParticleToHitListMap);

    pandora::StatusCode DivideCaloHits(const TwoDSlidingFitResult &slidingFitResult, const pandora::CartesianVector &splitPosition,
                                       pandora::CaloHitList &firstCaloHitList, pandora::CaloHitList &secondCaloHitList) const;

    pandora::StatusCode SplitCluster(const pandora::Cluster *const pCluster, const pandora::CaloHitList &caloHitList1,
                                     const pandora::CaloHitList &caloHitList2, pandora::ClusterVector &clusterSplittingList) const;

    std::string m_caloHitListName;
    std::string m_clusterListName;
    std::string m_mcParticleListName;
    std::string m_secVertexListName;
    float m_minFractionMerged;
    float m_minSecVertexAccuracy;
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_CLUSTER_SPLITTING_ALGORITHM_H
