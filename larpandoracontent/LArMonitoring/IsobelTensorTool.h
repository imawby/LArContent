/**
 *  @file   larpandoracontent/LArMonitoring/IsobelTensorTool.h
 *
 *  @brief  Header file for the isobel tensor tool class.
 *
 *  $Log: $
 */
#ifndef ISOBEL_TENSOR_TOOL_H
#define ISOBEL_TENSOR_TOOL_H 1

#include "larpandoracontent/LArThreeDReco/LArTransverseTrackMatching/ThreeDTransverseTracksAlgorithm.h"

namespace lar_content
{

/**
 *  @brief  IsobelTensorTool class
 */
class IsobelTensorTool : public TransverseTensorTool
{
public:
    /**
     *  @brief  Default constructor
     */
    IsobelTensorTool();

    bool Run(ThreeDTransverseTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor);

private:
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    float m_minXOverlapFraction;

};

} // namespace lar_content

#endif // #ifndef CLEAR_TRACKS_TOOL_H
