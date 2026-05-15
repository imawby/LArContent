/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerMatching/DLTwoViewMergeAndCreateShowersTool.h
 *
 *  @brief  Header file for the clear showers tool class.
 *
 *  $Log: $
 */
#ifndef DL_TWO_VIEW_MERGE_AND_CREATE_SHOWERS_TOOL_H
#define DL_TWO_VIEW_MERGE_AND_CREATE_SHOWERS_TOOL_H 1

#include "larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLMultiViewMatchingAlgorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  DLTwoViewMergeAndCreateShowersTool class
 */
class DLTwoViewMergeAndCreateShowersTool : public DLShowerMatchingTool
{
public:
    /**
     *  @brief  Default constructor
     */
    DLTwoViewMergeAndCreateShowersTool();

    bool Run(DLMultiViewMatchingAlgorithm *const pAlgorithm,
        const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


    bool CreateAmbiguousShower(DLMultiViewMatchingAlgorithm *const pAlgorithm,
        const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup, const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix,
        const pandora::HitType hitType1, const pandora::HitType hitType2);

    bool FindSeed(const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup,
        const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix, const pandora::HitType hitType1,
        const pandora::HitType hitType2, const pandora::Cluster *&pSeed1, const pandora::Cluster *&pSeed2);


    void MergeClusters(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup,
        const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix, const pandora::Cluster *const pSeed1,
        const pandora::Cluster *const pSeed2);

    float m_matchThreshold;   ///< Threshold score for match
};

} // namespace lar_dl_content

#endif // #ifndef DL_TWO_VIEW_MERGE_AND_CREATE_SHOWERS_TOOL_H
