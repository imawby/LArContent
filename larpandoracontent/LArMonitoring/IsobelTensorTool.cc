/**
 *  @file   larpandoracontent/LArMonitoring/IsobelTensorTool.cc
 *
 *  @brief  Implementation of the isobel tensor tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArMonitoring/IsobelTensorTool.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{

IsobelTensorTool::IsobelTensorTool() :
    m_minXOverlapFraction(5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool IsobelTensorTool::Run(ThreeDTransverseTracksAlgorithm *const pAlgorithm, TensorType &overlapTensor) {

    if (PandoraContentApi::GetSettings(*pAlgorithm)->ShouldDisplayAlgorithmInfo())
       std::cout << "----> Running Algorithm Tool: " << this->GetInstanceName() << ", " << this->GetType() << std::endl;

    std::cout << "Hello!" << std::endl;

    // Gets a list of significant U clusters
    ClusterVector sortedKeyClusters;
    overlapTensor.GetSortedKeyClusters(sortedKeyClusters);

    for(const Cluster * pCluster : sortedKeyClusters) 
    {
        // Get the elements for each of the U clusters
        TensorType::ElementList elementList;
        overlapTensor.GetConnectedElements(pCluster, true, elementList);



        std::cout << LArClusterHelper::GetClusterHitType(pCluster) << std::endl;
    }

    return false;

}

  
StatusCode IsobelTensorTool::ReadSettings(const TiXmlHandle xmlHandle)
{

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinXOverlapFraction", m_minXOverlapFraction));

    return STATUS_CODE_SUCCESS;
}
  

} // namespace lar_content
