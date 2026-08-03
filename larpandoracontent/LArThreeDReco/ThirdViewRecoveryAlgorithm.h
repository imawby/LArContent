/**
 *  @file   larpandoracontent/LArCheating/ThirdViewRecoveryAlgorithm.h
 *
 *  @brief  Header file for the challenge track shower if class.
 *
 *  $Log: $
 */
#ifndef LAR_THIRD_VIEW_RECOVERY_ALGORITHM_H
#define LAR_THIRD_VIEW_RECOVERY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

namespace lar_content
{

/**
 *  @brief  ThirdViewRecoveryAlgorithm class
 */
class ThirdViewRecoveryAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ThirdViewRecoveryAlgorithm();

private:
    typedef std::vector<pandora::HitType> HitTypes;
    
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void Proccess(const pandora::Pfo *const pPfo);

    pandora::StatusCode GetThirdViewProjection(const pandora::Cluster *const pCluster1, const pandora::Cluster *const pCluster2,
        float &minX, float &maxX, pandora::CartesianPointVector &projection);

    void GetThirdViewProjection(const pandora::CaloHitList &caloHitList, const TwoDSlidingFitResult &slidingFit1,
        const TwoDSlidingFitResult &slidingFit2, const float &xMin, const float &xMax, pandora::CartesianPointVector &projection);

    void GetMatchedCluster(const pandora::CartesianPointVector &projection, const std::vector<pandora::HitType> &hitTypes,
        const pandora::Cluster *&pClosestCluster);
    
    bool PassQualityCuts(const pandora::CartesianPointVector &projections, const pandora::Cluster *const pMatchedCluster);
    
    void ProcessAvailable(const pandora::Cluster *const pMatchedCluster, const pandora::CartesianPointVector &projection,
        const float minX, const float maxX, pandora::CaloHitList &collectedHits);

    void RecoverHitsWithoutMatchedClusterFit(const pandora::Cluster *const pMatchedCluster, const pandora::CartesianPointVector &projections,
        const float minX, const float maxX, pandora::CaloHitList &foundCaloHitList);
    
    void RecoverHitsWithMatchedClusterFit(const pandora::Cluster *const pMatchedCluster, const pandora::CartesianPointVector &matchedProjections,
        const pandora::CartesianPointVector &projections, const float minX, const float maxX, pandora::CaloHitList &foundCaloHitList);

    void ProcessThreeView(const pandora::Cluster *const pMatchedCluster, const pandora::CartesianPointVector &projection,
        const std::vector<pandora::HitType> &hitTypes, const float minX, const float maxX, pandora::CaloHitList &collectedHits);

    

    void GetProjectionInRange(const pandora::CaloHitList &caloHitList1, const pandora::CaloHitList &caloHitList2, pandora::CartesianPointVector &projections);
    



    void ProcessTwoView(const pandora::Cluster *const pMatchedCluster, const pandora::CartesianPointVector &projection,
        const std::vector<pandora::HitType> &hitTypes, const float minX, const float maxX, pandora::CaloHitList &collectedHits);    

    const pandora::Cluster* Get2DCluster(const pandora::Pfo *const pPfo, const pandora::HitType &hitType);



    void SplitIntoHitsAndIsolated(const pandora::Cluster *const pMatchedCluster, pandora::CaloHitList &collectedHits, pandora::CaloHitList &isolatedCollectedHits);

    void GetParentPfo(const pandora::Cluster *const pMatchedCluster, const pandora::Pfo *&pMatchedPfo);



    void ReassignHits(const pandora::Pfo *const pPfoToRecover, const pandora::CaloHitList &collectedHits,
        const pandora::Cluster *const pMatchedCluster);

    void ProcessRemnant(const pandora::Pfo *const pMatchedPfo, const pandora::Cluster *const pMatchedCluster);
    
    std::vector<pandora::HitType> GetViews(const pandora::Pfo *const pPfo);


    std::string m_trackPfoListName; ///< The name of the input track pfo list
    std::string m_showerPfoListName; ///< The name of the input shower pfo list
    
    std::string m_caloHitListNameU;   ///< The input calo hit list name for the U view
    std::string m_caloHitListNameV;   ///< The input calo hit list name for the V view
    std::string m_caloHitListNameW;   ///< The input calo hit list name for the W view    


    
    std::string m_clusterListNameU;   ///< The input cluster list name for the U view
    std::string m_clusterListNameV;   ///< The input cluster list name for the V view
    std::string m_clusterListNameW;   ///< The input cluster list name for the W view    
    unsigned int m_minNCaloHits;
    int m_slidingFitWindow;
    float m_matchedClusterMaxSep;
    float m_gapTolerance;
    float m_minMatchedFrac;
    float m_recoveryMaxTransSep;
    float m_keepMaxTransSep;
    float m_matchedXRange;
    
};

} // namespace lar_content

#endif // #ifndef LAR_THIRD_VIEW_RECOVERY_ALGORITHM_H
