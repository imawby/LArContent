/**
 *  @file   larpandoracontent/LArThreeDReco/LArShowerMatching/DLTwoViewClearShowersTool.h
 *
 *  @brief  Header file for the clear showers tool class.
 *
 *  $Log: $
 */
#ifndef DL_TWO_VIEW_CLEAR_SHOWERS_TOOL_H
#define DL_TWO_VIEW_CLEAR_SHOWERS_TOOL_H 1

#include "larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLMultiViewMatchingAlgorithm.h"

namespace lar_dl_content
{

/**
 *  @brief  DLTwoViewClearShowersTool class
 */
class DLTwoViewClearShowersTool : public DLShowerMatchingTool
{
public:
    /**
     *  @brief  Default constructor
     */
    DLTwoViewClearShowersTool();

    bool Run(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::SimilarityMatrix &globalSimMatrix);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


    void CreateClearShowers(DLMultiViewMatchingAlgorithm *const pAlgorithm, const DLMultiViewMatchingAlgorithm::ClusterGroup &clusterGroup);
};

} // namespace lar_dl_content

#endif // #ifndef DL_THREE_VIEW_CLEAR_SHOWERS_TOOL_H
