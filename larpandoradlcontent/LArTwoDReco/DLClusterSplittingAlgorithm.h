/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/CleanDLClusterSplittingAlgorithm.h
 *
 *  @brief  Header file for the two dimensional sliding fit splitting algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_CLEAN_DL_CLUSTER_SPLITTING_ALGORITHM_H
#define LAR_CLEAN_DL_CLUSTER_SPLITTING_ALGORITHM_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArUtility/KalmanFilter.h"
#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"
#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/ClusterSplittingAlgorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  CleanDLClusterSplittingAlgorithm class
 */
class CleanDLClusterSplittingAlgorithm : public pandora::Algorithm
{

class ClusterHit
{
public:

    ClusterHit(const pandora::CaloHit *const pCaloHit, const float l, const float t);

    /**
     *  @brief  HierarchyPfo == operator
     *
     *  @param  rhs the pfo to compare
     */
    bool operator==(const ClusterHit &rhs) const;

    float m_l;
    float m_t;
    const pandora::CaloHit *m_pHit;
};

class Feature
{
public:

    Feature(const float mean, const float std, const pandora::FloatVector &sequence);
    void Normalise();
    void Smooth();

    float m_mean;
    float m_std;
    pandora::FloatVector m_sequence;
};

public:
    /**
     *  @brief  Default constructor
     */
    CleanDLClusterSplittingAlgorithm();

    pandora::StatusCode Run();

private:
    typedef std::vector<ClusterHit> ClusterPath;
    typedef lar_content::KDTreeLinkerAlgo<const pandora::CaloHit *, 2> HitKDTree2D;
    typedef lar_content::KDTreeNodeInfoT<const pandora::CaloHit *, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;
    //typedef std::map<int, std::pair<const pandora::CaloHit*, float>> ClusterPath;
    typedef std::map<std::string, Feature> Features;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    pandora::StatusCode GetLists();

    void ProcessCluster(const pandora::Cluster *const pCluster, HitKDTree2D &kdTree);

    void FindPath(const pandora::Cluster *const pCluster, const lar_content::TwoDSlidingFitResult &clusterFit, ClusterPath &clusterPath) const;

    void InitialiseFeatures(Features &features) const;

    void FillFeatures(const ClusterPath &clusterPath, const lar_content::TwoDSlidingFitResult &clusterFit, Features &features, HitKDTree2D &kdTree) const;

    lar_content::KalmanFilter2D InitialiseKalmanFilter(const ClusterPath &clusterPath) const;

    float GetAngle(const pandora::CartesianVector &position2D, const lar_content::TwoDSlidingFitResult &clusterFit) const;

    float GetDistanceToGap(const pandora::CaloHit *const pPrevHit, const pandora::CaloHit *const pCurrentHit) const;

    float GetDistanceToEventHit(HitKDTree2D &kdTree, const pandora::CaloHit *const pCaloHit, const pandora::CaloHitList &clusterPathHits) const;

    float GetDistanceToClusterHit(const ClusterHit &currentHit, const ClusterHit &nextHit) const;

    float GetDistanceToSecVertex(const pandora::CaloHit *const pCaloHit, const pandora::HitType hitType) const;

    void GetSplitIndices(const Features &features, pandora::IntVector &splitIndices);

    void GetSplitIndices(const Features &features, const pandora::IntVector &windowStart, 
        pandora::IntVector &splitIndices, pandora::FloatVector &splitScores);

    void FilterModelOutput(const pandora::FloatVector &splitScores, pandora::IntVector &splitIndices) const;

    void FilterSplitIndices(const ClusterPath &clusterPath, const pandora::HitType hitType, pandora::IntVector &splitIndices) const;

    std::vector<pandora::CaloHitList> DivideCaloHits(const pandora::CaloHitList &clusterHits, const ClusterPath &clusterPath, 
        const lar_content::TwoDSlidingFitResult &clusterFit, const pandora::IntVector &splitIndices) const;

    void SplitCluster(const pandora::Cluster *const pCluster, const std::vector<pandora::CaloHitList> &splitClusterHits) const;

    std::string m_caloHitListName;
    std::string m_clusterListName;
    std::string m_nuVertexListName;
    std::string m_secVertexListName;
    const pandora::CaloHitList *m_pCaloHitList;
    const pandora::ClusterList *m_pClusterList;
    const pandora::VertexList *m_pNuVertexList;
    const pandora::VertexList *m_pSecVertexList;
    unsigned int m_minClusterHits;
    int m_slidingWindow;
    float m_lBinSize;
    float m_searchRegion1D;
    unsigned int m_windowLength;
    float m_isContaminatedThreshold;
    float m_isSplitThreshold;

    std::string m_windowModelName;
    LArDLHelper::TorchModel m_windowModel;
    std::string m_splitPosModelName;
    LArDLHelper::TorchModel m_splitPosModel;
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline CleanDLClusterSplittingAlgorithm::ClusterHit::ClusterHit(const pandora::CaloHit *const pCaloHit, const float l, const float t) :
    m_l(l),
    m_t(t),
    m_pHit(pCaloHit)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline bool CleanDLClusterSplittingAlgorithm::ClusterHit::operator==(const ClusterHit &rhs) const
{
    return this->m_l == rhs.m_l;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline CleanDLClusterSplittingAlgorithm::Feature::Feature(const float mean, const float std, const pandora::FloatVector &sequence) :
    m_mean(mean),
    m_std(std),
    m_sequence(sequence)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void CleanDLClusterSplittingAlgorithm::Feature::Normalise()
{
    for (float &val : m_sequence)
        val = (val - m_mean) / m_std;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void CleanDLClusterSplittingAlgorithm::Feature::Smooth()
{
    pandora::FloatVector temp(m_sequence);
    
    for (unsigned int iEntry = 0; iEntry < m_sequence.size(); ++iEntry)
    {
        float total(temp.at(iEntry));
        int nEntries(1);

        for (unsigned int iWindow = 1; iWindow <= 4; ++iWindow)
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

        m_sequence[iEntry] = (total / nEntries);
    }
}

} // namespace lar_content

#endif // #ifndef LAR_CLEAN_DL_CLUSTER_SPLITTING_ALGORITHM_H
