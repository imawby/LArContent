/**
 *  @file   larpandoracontent/LArMonitoring/VertexDrawingAlgorithm.cc
 *
 *  @brief  Implementation of the vertex drawing algorithm
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArMonitoring/VertexDrawingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

VertexDrawingAlgorithm::VertexDrawingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexDrawingAlgorithm::Run()
{
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

    const VertexList *pNuVertexList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_nuVertexListName, pNuVertexList));

    if (!pNuVertexList || (pNuVertexList->size() != 1))
    {
        std::cout << "NO NEUTRINO VERTEX WAS RECONSTRUCTED, NOTHING TO VISUALIZE!" << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    const CartesianVector &nuVertex3D = pNuVertexList->front()->GetPosition();
    const CartesianVector nuVertexU = LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, TPC_VIEW_U);
    const CartesianVector nuVertexV = LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, TPC_VIEW_V);
    const CartesianVector nuVertexW = LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, TPC_VIEW_W);

    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &nuVertexU, "U reco nu vertex", BLUE, 2));
    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &nuVertexV, "V reco nu vertex", RED, 2));
    PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &nuVertexW, "W reco nu vertex", BLACK, 2));

    PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode VertexDrawingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NuVertexListName", m_nuVertexListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
