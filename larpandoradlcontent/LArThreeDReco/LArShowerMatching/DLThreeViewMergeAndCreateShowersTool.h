/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerMatching/DLThreeViewMergeAndCreateShowersTool.h
 *
 *  @brief  Header file for the clear showers tool class.
 *
 *  $Log: $
 */
#ifndef DL_THREE_VIEW_MERGE_AND_CREATE_SHOWERS_TOOL_H
#define DL_THREE_VIEW_MERGE_AND_CREATE_SHOWERS_TOOL_H 1

#include "larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLMultiViewMatchingAlgorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  DLThreeViewMergeAndCreateShowersTool class
 */
class DLThreeViewMergeAndCreateShowersTool : public DLShowerMatchingTool
{
public:
    /**
     *  @brief  Default constructor
     */
    DLThreeViewMergeAndCreateShowersTool();

    bool Run(DLMultiViewMatchingAlgorithm *const pAlgorithm,
        const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


    bool CreateAmbiguousShower(DLMultiViewMatchingAlgorithm *const pAlgorithm,
        const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup, const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix);

    bool FindSeed(const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup,
        const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix, const pandora::Cluster *&pSeedU,
        const pandora::Cluster *&pSeedV, const pandora::Cluster *&pSeedW);


    void MergeClusters(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup,
        const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix, const pandora::Cluster *const pSeedU,
        const pandora::Cluster *const pSeedV, const pandora::Cluster *const pSeedW);

    float m_matchThreshold;   ///< Threshold score for match
};

} // namespace lar_dl_content

#endif // #ifndef DL_THREE_VIEW_MERGE_AND_CREATE_SHOWERS_TOOL_H
