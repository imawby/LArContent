/**
 *  @file   larpandoracontent/LArCheating/CheatingKalmanSplittingAlgorithm.h
 *
 *  @brief  Header file for the cheating cluster splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_KALMAN_SPLITTING_ALGORITHM_H
#define LAR_CHEATING_KALMAN_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  CheatingKalmanSplittingAlgorithm class
 */
class CheatingKalmanSplittingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingKalmanSplittingAlgorithm();

    ~CheatingKalmanSplittingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    typedef std::unordered_map<const pandora::MCParticle *, pandora::CaloHitList> MCParticleToHitListMap;
    typedef std::unordered_map<const pandora::CaloHit *, const pandora::MCParticle *> HitToMCParticleMap;
    typedef std::unordered_map<const pandora::Cluster *, const pandora::MCParticle *> ClusterToMCParticleMap;
    typedef std::unordered_map<const pandora::Cluster *, std::vector<std::pair<const pandora::MCParticle*, pandora::CaloHitList>>> ClusterToMCParticleListMap;
    typedef std::map<int, std::pair<const pandora::CaloHit*, float>> ClusterPath;
    typedef KDTreeLinkerAlgo<const pandora::CaloHit *, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit *, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;

    void FillPandoraMaps(const pandora::ClusterList *const pClusterList, const pandora::CaloHitList *const pCaloHitList, 
                         const pandora::MCParticleList *const pMCParticleList, const pandora::VertexList *const pSecVertexList,
                         MCParticleToHitListMap &mcParticleToHitListMap, HitToMCParticleMap &hitToMCParticleMap, ClusterToMCParticleMap &clusterToMCParticleMap,
                         ClusterToMCParticleListMap &clusterToMCParticleListMap);

    void ProbeContaminants(const pandora::ClusterList *const pClusterList, const pandora::CaloHitList *const pCaloHitList, 
                           const pandora::VertexList *const pSecVertexList, ClusterToMCParticleMap &clusterToMCParticleMap, 
                           ClusterToMCParticleListMap &clusterToMCParticleListMap, MCParticleToHitListMap &mcParticleToHitListMap);

    void FindPath(const pandora::Cluster *const pCluster, const TwoDSlidingFitResult &clusterFit, ClusterPath &clusterPath);

    void PerformKalmanFit(const ClusterPath &clusterPath, const float totalEnergy, std::vector<float> &posDiffDist, std::vector<float> &energyDiffDist, 
                          std::vector<float> &scatterDist, std::vector<float> &mahalanobisDist);

    void BuildKDTree(const pandora::Cluster *const pCluster, const pandora::CaloHitList *const pCaloHitList, HitKDTree2D &kdTree_event);

    float GetDistanceToEventHit(HitKDTree2D &kdTree_event, const pandora::Cluster *const pCluster, 
                                const pandora::CaloHit *const pCaloHit);

    float GetDistanceToClusterHit(const ClusterPath::iterator &currentHit, const ClusterPath::iterator &endIter);

    float GetDistanceToGap(const pandora::CartesianVector &position2D) const;

    float GetDistanceToSecVertex(const pandora::CaloHit *const pCaloHit, const pandora::VertexList *const pSecVertexList, 
                                 const pandora::HitType hitType);


    std::string m_caloHitListName;
    std::string m_clusterListName;
    std::string m_mcParticleListName;
    std::string m_secVertexListName;
    unsigned int m_minClusterHits;
    int m_minTargetMCHits;
    float m_minFractionMerged;
    int m_slidingWindow;
    float m_lBinSize;
    float m_endpointBuffer;
    float m_searchRegion1D;
    bool m_writeVisInfo;
    std::string m_treeName;
    std::string m_fileName;
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_KALMAN_SPLITTING_ALGORITHM_H
