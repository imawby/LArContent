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

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

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
    typedef std::unordered_map<const pandora::Cluster *, pandora::MCParticleList> ClusterToMCParticleListMap;
    typedef std::unordered_map<const pandora::MCParticle *, pandora::CartesianVector> MCParticleSecVertexMap;

    void FillPandoraMaps(const pandora::ClusterList *const pClusterList, const pandora::CaloHitList *const pCaloHitList, 
                         const pandora::MCParticleList *const pMCParticleList, const pandora::VertexList *const pSecVertexList,
                         MCParticleToHitListMap &mcParticleToHitListMap, HitToMCParticleMap &hitToMCParticleMap, ClusterToMCParticleMap &clusterToMCParticleMap,
                         ClusterToMCParticleListMap &clusterToMCParticleListMap, MCParticleSecVertexMap &mcParticleSecVertexMap);

    void ProbeContaminants(const pandora::ClusterList *const pClusterList, const pandora::CaloHitList *const pCaloHitList, 
                           const pandora::VertexList *const pSecVertexList, ClusterToMCParticleMap &clusterToMCParticleMap, 
                           ClusterToMCParticleListMap &clusterToMCParticleListMap, MCParticleToHitListMap &mcParticleToHitListMap);

    void FindPath(const pandora::Cluster *const pCluster, const TwoDSlidingFitResult &clusterFit, 
                  std::map<int, std::pair<const pandora::CaloHit*, float>> &clusterPath);

    void PerformKalmanFit(const std::map<int, std::pair<const pandora::CaloHit*, float>> &clusterPath, 
                          std::vector<float> &posDiffDist, std::vector<float> &energyDiffDist, std::vector<float> &widthDiffDist, 
                          std::vector<float> &scatterDist, std::vector<float> &mahalanobisDist);

    float GetDistanceToGap(const pandora::CartesianVector &position2D) const;

    float GetDistanceToEventHit(const pandora::CaloHit *const pCaloHit, const pandora::CaloHitList &clusterHits, 
                                const pandora::CaloHitList *const pEventHits);

    std::string m_caloHitListName;
    std::string m_clusterListName;
    std::string m_mcParticleListName;
    std::string m_secVertexListName;
    float m_minFractionMerged;
    float m_minSecVertexAccuracy;
    bool m_writeFile;
    std::string m_treeName;
    std::string m_fileName;
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_KALMAN_SPLITTING_ALGORITHM_H
