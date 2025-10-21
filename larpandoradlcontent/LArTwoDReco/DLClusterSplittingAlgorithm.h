/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/DLClusterSplittingAlgorithm.h
 *
 *  @brief  Header file for the two dimensional sliding fit splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_CLUSTER_SPLITTING_ALGORITHM_H
#define LAR_DL_CLUSTER_SPLITTING_ALGORITHM_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"
#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/ClusterSplittingAlgorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  DLClusterSplittingAlgorithm class
 */
class DLClusterSplittingAlgorithm : public pandora::Algorithm
{

class Features
{
public:

    float NormaliseTransverse(const float feature);
    float NormaliseEnergy(const float feature);
    float NormaliseHitWidth(const float feature);
    float NormaliseTheta(const float feature);
    float NormaliseSecVertex(const float feature);
    float NormaliseGapSep(const float feature);
    float NormaliseEventHitSep(const float feature);
    float NormaliseClusterHitSep(const float feature);
    void SmoothFeature(pandora::FloatVector &feature);

    pandora::FloatVector m_transverse;
    pandora::FloatVector m_energy;
    pandora::FloatVector m_hitWidth;
    pandora::FloatVector m_theta;
    pandora::FloatVector m_secVertex;
    pandora::FloatVector m_gapSep;
    pandora::FloatVector m_eventHitSep;
    pandora::FloatVector m_clusterHitSep;
};

public:
    /**
     *  @brief  Default constructor
     */
    DLClusterSplittingAlgorithm();

    pandora::StatusCode Run();

private:
    typedef lar_content::KDTreeLinkerAlgo<const pandora::CaloHit *, 2> HitKDTree2D;
    typedef lar_content::KDTreeNodeInfoT<const pandora::CaloHit *, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;
    typedef std::map<int, std::pair<const pandora::CaloHit*, float>> ClusterPath;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StatusCode ThisDivideCaloHits(const pandora::Cluster *const pCluster);

    void FindPath(const pandora::Cluster *const pCluster, const lar_content::TwoDSlidingFitResult &clusterFit, ClusterPath &clusterPath) const;

    void FillFeatures(const ClusterPath &clusterPath, const lar_content::TwoDSlidingFitResult &clusterFit, Features &features) const;

    float GetAngle(const pandora::CartesianVector &position2D, const lar_content::TwoDSlidingFitResult &clusterFit) const;

    float GetDistanceToGap(const pandora::CartesianVector &position2D, const pandora::HitType hitType) const;

    float GetDistanceToEventHit(HitKDTree2D &kdTree, const pandora::CaloHit *const pCaloHit, const pandora::CaloHitList &clusterPathHits) const;

    float GetDistanceToClusterHit(const ClusterPath::const_iterator &currentHit, const ClusterPath::const_iterator &endIter) const;

    float GetDistanceToSecVertex(const pandora::CaloHit *const pCaloHit, const pandora::HitType hitType) const;

    void GetWindows(Features &features, pandora::IntVector &splitIndices);

    std::string m_caloHitListName;
    std::string m_clusterListName;
    std::string m_secVertexListName;
    const pandora::VertexList *m_pSecVertexList;
    unsigned int m_minClusterHits;
    int m_slidingWindow;
    float m_lBinSize;
    float m_searchRegion1D;
    std::string m_windowModelName;
    LArDLHelper::TorchModel m_windowModel;
    std::string m_splitPosModelName;
    LArDLHelper::TorchModel m_splitPosModel;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DLClusterSplittingAlgorithm::Features::NormaliseTransverse(const float feature)
{
    /* mean: -0.007181927387650833 */
    /* stan_dev: 2.5722720356088034 */
    return (feature - (-0.01)) / 2.57;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DLClusterSplittingAlgorithm::Features::NormaliseEnergy(const float feature)
{
    /* mean: 0.48753944268089333 */
    /* stan_dev: 0.2795828223860681 */
    return (feature - 0.49) / 0.28;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DLClusterSplittingAlgorithm::Features::NormaliseHitWidth(const float feature)
{
    /* mean: 0.5638581256164427 */
    /* stan_dev: 0.16363261394475448 */
    return (feature - 0.56) / 0.16;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DLClusterSplittingAlgorithm::Features::NormaliseTheta(const float feature)
{
    /* mean: -0.011612029556346253 */
    /* stan_dev: 0.10516526130812948 */
    return (feature - (-0.01)) / 0.11;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DLClusterSplittingAlgorithm::Features::NormaliseSecVertex(const float feature)
{
    /* mean: 93.18852645417384 */
    /* stan_dev: 103.39096393057723 */
    return (feature - 93.19) / 103.39;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DLClusterSplittingAlgorithm::Features::NormaliseGapSep(const float feature)
{
    /* mean: 171.77832424311225 */
    /* stan_dev: 101.44208765041625 */
    return (feature - 171.78) / 101.44;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DLClusterSplittingAlgorithm::Features::NormaliseEventHitSep(const float feature)
{
    /* mean: 5.003103494406862 */
    /* stan_dev: 5.832603908694449 */
    return (feature - 5.00) / 5.83;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float DLClusterSplittingAlgorithm::Features::NormaliseClusterHitSep(const float feature)
{
    /* mean: 0.5379566560251992 */
    /* stan_dev: 0.1261325429847093 */
    return (feature - 0.54) / 0.13;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void DLClusterSplittingAlgorithm::Features::SmoothFeature(pandora::FloatVector &feature)
{
    pandora::FloatVector temp(feature);
    
    for (unsigned int iEntry = 0; iEntry < feature.size(); ++iEntry)
    {
        float total(0.f);
        int nEntries(0);

        for (unsigned int iWindow = 0; iWindow < 5; ++iWindow)
        {
            const int belowIndex(iEntry - iWindow);
            const unsigned int aboveIndex(iEntry + iWindow);
            
            if (belowIndex >= 0)
            {
                total += temp.at(belowIndex);
                ++nEntries;
            }
                
            if (aboveIndex < temp.size())
            {
                total += temp.at(aboveIndex);
                ++nEntries;
            }
        }

        feature[iEntry] = (total / nEntries);
    }
}

} // namespace lar_content

#endif // #ifndef LAR_DL_CLUSTER_SPLITTING_ALGORITHM_H
