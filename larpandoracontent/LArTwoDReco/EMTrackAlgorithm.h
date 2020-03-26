/**
 *  @file   EMTrackAlgorithm.h
 *
 *  @brief  Header file for the em track algorithm class
 *
 *  $Log: $
 */
#ifndef LAR_EM_TRACK_ALGORITHM_H
#define LAR_EM_TRACK_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{


/**
 *  @brief  EMTrackAlgorithm class
 */
class EMTrackAlgorithm : public pandora::Algorithm
{
public:

    class ClusterAssociation
    {
    public:
        ClusterAssociation(const pandora::Cluster *const pAssociatedCluster, const pandora::CartesianVector &innerMergePoint, const pandora::CartesianVector &innerMergeDirection, const pandora::CartesianVector &outerMergePoint, const pandora::CartesianVector &outerMergeDirection);

        const pandora::Cluster* GetAssociatedCluster() const;
        pandora::CartesianVector GetInnerMergePoint() const;
        pandora::CartesianVector GetInnerMergeDirection() const;
        pandora::CartesianVector GetOuterMergePoint() const;
        pandora::CartesianVector GetOuterMergeDirection() const;
    private:
        const pandora::Cluster*     m_pAssociatedCluster;
        pandora::CartesianVector    m_innerMergePoint;
        pandora::CartesianVector    m_innerMergeDirection;
        pandora::CartesianVector    m_outerMergePoint;
        pandora::CartesianVector    m_outerMergeDirection;
    };

    /**
     *  @brief  Default constructor
     */
    EMTrackAlgorithm();

    typedef std::unordered_map<const pandora::Cluster*, const pandora::Cluster*> ClusterToAssociatedClusterMap;
    typedef std::unordered_map<const pandora::CaloHit*, const pandora::Cluster*> CaloHitToParentClusterMap;
    typedef std::unordered_map<const pandora::Cluster*, ClusterAssociation> ClusterToClusterAssociationMap;
    
 private:
    pandora::StatusCode Run();
    
    void SelectCleanClusters(const pandora::ClusterList *pClusterList, pandora::ClusterVector &clusterVector);

    void InitialiseSlidingFitResultMap(const pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitResultMap);

    void UpdateSlidingFitResultMap(const pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitResultMap);

    void MergeClusters(const pandora::Cluster *const pClusterToEnlarge, const pandora::Cluster *const pClusterToDelete, pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitResultMap);

    bool FindBestClusterAssociation(const pandora::Cluster *const pCurrentCluster, const pandora::ClusterVector &clusterVector, const TwoDSlidingFitResultMap &slidingFitResultMap, ClusterAssociation &clusterAssociation);
    
    bool AreClustersAssociated(const pandora::CartesianVector &currentPoint, const pandora::CartesianVector &currentDirection, const pandora::CartesianVector &testPoint, const pandora::CartesianVector &testDirection);
    
    void GetClusterMergingCoordinates(const TwoDSlidingFitResult &cluster1FitResult, pandora::CartesianVector &cluster1MergePoint, pandora::CartesianVector &cluster1Direction, const TwoDSlidingFitResult &cluster2FitResult, pandora::CartesianVector &cluster2MergePoint, pandora::CartesianVector &cluster2Direction);
    
    void GetExtrapolatedCaloHits(const pandora::Cluster *const innerCluster, const ClusterAssociation &clusterAssociation, const pandora::CaloHitList *const pCaloHitList, const pandora::ClusterList *const pClusterList, pandora::CaloHitVector &extrapolatedCaloHitVector, CaloHitToParentClusterMap &caloHitToParentClusterMap);

    bool IsTrackContinuous(const ClusterAssociation &clusterAssociation, pandora::CaloHitVector &extrapolatedCaloHitVector);
    
    void ConnectByLine(const pandora::CartesianVector &innerCoordinate, const pandora::CartesianVector &outerCoordinate, float &gradient, float &zIntercept);

    void AddHitsToCluster(const pandora::Cluster *const clusterToEnlarge, pandora::ClusterVector &clusterVector, TwoDSlidingFitResultMap &slidingFitResultMap, const CaloHitToParentClusterMap &caloHitToParentClusterMap, const pandora::CaloHitVector &extrapolatedCaloHitVector);
    
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;
    CaloHitToParentClusterMap m_caloHitToParentClusterMap;
    
    float m_minCaloHits; // minimum number of calo hits for cluster to be considered
    float m_maxXSeparation;
    float m_maxZSeparation;

};

//------------------------------------------------------------------------------------------------------------------------------------------

inline const pandora::Cluster* EMTrackAlgorithm::ClusterAssociation::GetAssociatedCluster() const
{
    return m_pAssociatedCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector EMTrackAlgorithm::ClusterAssociation::GetInnerMergePoint() const
{
    return m_innerMergePoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector EMTrackAlgorithm::ClusterAssociation::GetInnerMergeDirection() const
{
    return m_innerMergeDirection;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector EMTrackAlgorithm::ClusterAssociation::GetOuterMergePoint() const
{
    return m_outerMergePoint;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline pandora::CartesianVector EMTrackAlgorithm::ClusterAssociation::GetOuterMergeDirection() const
{
    return m_outerMergeDirection;
}

} //namespace lar_content


#endif // #ifndef LAR_EM_TRACK_ALGORITHM_H