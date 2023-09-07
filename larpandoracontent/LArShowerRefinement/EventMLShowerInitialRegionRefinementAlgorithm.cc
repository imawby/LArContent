/**
 *  @file   larpandoracontent/LArShowerRefinement/EventMLShowerInitialRegionRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the electron initial region refinement algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingEventShowerStartFinderTool.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"
#include "larpandoracontent/LArShowerRefinement/EventMLShowerInitialRegionRefinementAlgorithm.h"

#ifdef MONITORING
#include "PandoraMonitoringApi.h"
#endif

using namespace pandora;

namespace lar_content
{

EventMLShowerInitialRegionRefinementAlgorithm::EventMLShowerInitialRegionRefinementAlgorithm() :
    m_minShowerHits3D(50)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventMLShowerInitialRegionRefinementAlgorithm::Run()
{
    //////////////////////////////////////////
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
    //////////////////////////////////////////

    // Get 2D shower starts
    CartesianPointVector showerStartsU, showerStartsV, showerStartsW;
    this->GetShowerStarts(showerStartsU, showerStartsV, showerStartsW);

    if (showerStartsU.empty() || showerStartsV.empty() || showerStartsW.empty())
        return STATUS_CODE_SUCCESS;

    // Get target reco showers
    PfoVector showerPfoVector;
    this->FillShowerPfoVector(showerPfoVector);

    if (showerPfoVector.empty())
        return STATUS_CODE_SUCCESS;

    for (const ParticleFlowObject *const pShowerPfo : showerPfoVector)
        this->RefineShower(pShowerPfo);

    PfoList visualisationPfos;
    visualisationPfos.insert(visualisationPfos.begin(), showerPfoVector.begin(), showerPfoVector.end());

    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &visualisationPfos, "ShowerPfos", RED);

    for (const CartesianVector &showerStart : showerStartsU)
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &showerStart, "showerStartU", BLACK, 2);

    for (const CartesianVector &showerStart : showerStartsV)
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &showerStart, "showerStartV", BLACK, 2);

    for (const CartesianVector &showerStart : showerStartsW)
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &showerStart, "showerStartW", BLACK, 2);

    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventMLShowerInitialRegionRefinementAlgorithm::GetShowerStarts(CartesianPointVector &showerStartsU, CartesianPointVector &showerStartsV, 
    CartesianPointVector &showerStartsW) const
{
    const MCParticleList *pMCParticleList(nullptr);

    if (PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList) != STATUS_CODE_SUCCESS)
        return;

    if (!pMCParticleList || pMCParticleList->empty())
    {
        std::cout << "EventMLShowerInitialRegionRefinementAlgorithm: unable to find mc particle list " << m_mcParticleListName << std::endl;
        return;
    }

    const CaloHitList *pCaloHitList(nullptr);

    if (PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList) != STATUS_CODE_SUCCESS)
        return;

    if (!pCaloHitList || pCaloHitList->empty())
    {
        std::cout << "EventMLShowerInitialRegionRefinementAlgorithm: unable to find calo hit list " << m_caloHitListName << std::endl;
        return;
    }

    m_pCheatingEventShowerStartFinderTool->Run(pMCParticleList, pCaloHitList, showerStartsU, showerStartsV, showerStartsW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventMLShowerInitialRegionRefinementAlgorithm::FillShowerPfoVector(PfoVector &showerPfoVector) const
{
    const PfoList *pPfoList(nullptr);

    if (PandoraContentApi::GetList(*this, m_showerPfoListName, pPfoList) != STATUS_CODE_SUCCESS)
        return;

    if (!pPfoList || pPfoList->empty())
    {
        std::cout << "EventMLShowerInitialRegionRefinementAlgorithm: unable to find shower pfo list " << m_showerPfoListName << std::endl;
        return;
    }

    for (const ParticleFlowObject *const pShowerPfo : *pPfoList)
    {
        // Only consider significant showers
        CaloHitList caloHits3D;
        LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHits3D);

        if (caloHits3D.size() < m_minShowerHits3D)
            continue;

        showerPfoVector.push_back(pShowerPfo);
    }

    std::sort(showerPfoVector.begin(), showerPfoVector.end(), LArPfoHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void EventMLShowerInitialRegionRefinementAlgorithm::RefineShower(const ParticleFlowObject *const pShowerPfo) const
{
    CartesianVector nuVertex3D(0.f, 0.f, 0.f);

    if (this->GetNeutrinoVertex(nuVertex3D) != STATUS_CODE_SUCCESS)
        return;



    /*
    // Create the 2D connetion pathways
    ProtoShowerVector protoShowerVectorU, protoShowerVectorV, protoShowerVectorW;

    this->BuildViewProtoShowers(pShowerPfo, nuVertex3D, TPC_VIEW_U, protoShowerVectorU);
    this->BuildViewProtoShowers(pShowerPfo, nuVertex3D, TPC_VIEW_V, protoShowerVectorV);
    this->BuildViewProtoShowers(pShowerPfo, nuVertex3D, TPC_VIEW_W, protoShowerVectorW);
    */
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventMLShowerInitialRegionRefinementAlgorithm::GetNeutrinoVertex(CartesianVector &nuVertex3D) const
{
    const VertexList *pNuVertexList(nullptr);
    const StatusCode statusCode(PandoraContentApi::GetList(*this, m_neutrinoVertexListName, pNuVertexList));

    if (statusCode != STATUS_CODE_SUCCESS)
        return statusCode;

    if (!pNuVertexList || (pNuVertexList->size() != 1))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "EventMLShowerInitialRegionRefinementAlgorithm: unable to find vertex list " << m_neutrinoVertexListName
                      << " if it does exist, it may have more than one nu vertex" << std::endl;

        return STATUS_CODE_NOT_INITIALIZED;
    }

    nuVertex3D = pNuVertexList->front()->GetPosition();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode EventMLShowerInitialRegionRefinementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoVertexListName", m_neutrinoVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinShowerHits3D", m_minShowerHits3D));

    AlgorithmTool *pAlgorithmTool(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "CheatingEventShowerStartFinder", pAlgorithmTool));
    m_pCheatingEventShowerStartFinderTool = dynamic_cast<CheatingEventShowerStartFinderTool *>(pAlgorithmTool);



    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
