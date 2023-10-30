/**
 *  @file   larpandoracontent/LArShowerRefinement/MLShowerInitialRegionRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the electron initial region refinement algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingShowerStartFinderTool.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"
#include "larpandoracontent/LArShowerRefinement/MLShowerInitialRegionRefinementAlgorithm.h"
#include "larpandoracontent/LArShowerRefinement/PeakDirectionFinderTool.h"
#include "larpandoracontent/LArShowerRefinement/ShowerSpineFinderTool.h"
#include "larpandoracontent/LArShowerRefinement/ShowerStartFinderTool.h"



#ifdef MONITORING
#include "PandoraMonitoringApi.h"
#endif

using namespace pandora;

namespace lar_content
{

MLShowerInitialRegionRefinementAlgorithm::MLShowerInitialRegionRefinementAlgorithm() :
    m_minShowerHits3D(50)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLShowerInitialRegionRefinementAlgorithm::Run()
{
    //////////////////////////////////////////
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
    //////////////////////////////////////////

    // Get target reco showers
    PfoVector showerPfoVector;
    this->FillShowerPfoVector(showerPfoVector);

    if (showerPfoVector.empty())
        return STATUS_CODE_SUCCESS;

    for (const ParticleFlowObject *const pShowerPfo : showerPfoVector)
        this->RefineShower(pShowerPfo);

    /*
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
    */

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLShowerInitialRegionRefinementAlgorithm::RefineShower(const ParticleFlowObject *const pShowerPfo) const
{
    CartesianVector nuVertex3D(0.f, 0.f, 0.f);

    if (this->GetNeutrinoVertex(nuVertex3D) != STATUS_CODE_SUCCESS)
        return;

    // Create the 2D connetion pathways
    ProtoShowerVector protoShowerVectorU, protoShowerVectorV, protoShowerVectorW;

    this->BuildViewProtoShowers(pShowerPfo, nuVertex3D, TPC_VIEW_U, protoShowerVectorU);
    this->BuildViewProtoShowers(pShowerPfo, nuVertex3D, TPC_VIEW_V, protoShowerVectorV);
    this->BuildViewProtoShowers(pShowerPfo, nuVertex3D, TPC_VIEW_W, protoShowerVectorW);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MLShowerInitialRegionRefinementAlgorithm::BuildViewProtoShowers(const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertex3D, HitType hitType, ProtoShowerVector &protoShowerVector) const
{
    const CaloHitList *pViewHitList(nullptr);

    if (this->GetHitListOfType(hitType, pViewHitList) != STATUS_CODE_SUCCESS)
        return;

    /*
    CartesianVector showerVertexPosition(0.f, 0.f, 0.f);
    try
    {
        showerVertexPosition = this->GetShowerVertex(pShowerPfo, hitType, nuVertex3D);
    }
    catch (...)
    {
        return;
    }
    */

    // Determine directions of pathways out of neutrino vertex
    CartesianPointVector peakDirectionVector;
    if (m_pShowerPeakDirectionFinderTool->Run(pShowerPfo, nuVertex3D, pViewHitList, hitType, peakDirectionVector) != STATUS_CODE_SUCCESS)
        return;

    // Investigate each direction
    CaloHitList unavailableHitList;
    for (CartesianVector &peakDirection : peakDirectionVector)
    {
        std::cout << "Finding associated shower spine..." << std::endl;
        // Collect the hits associated with the pathway (the shower spine)
        CaloHitList showerSpineHitList;
        if (m_pShowerSpineFinderTool->Run(nuVertex3D, pViewHitList, hitType, peakDirection, unavailableHitList, showerSpineHitList) != STATUS_CODE_SUCCESS)
        {
            std::cout << "No spine found!" << std::endl;
            continue;
        }

        /*
        this->RefineShowerVertex(pShowerPfo, hitType, nuVertex3D, peakDirection, showerVertexPosition);

        // Quality check: If the spine passes the shower vertex, does it live inside the shower?
        if (!this->IsSpineCoincident(pShowerPfo, nuVertex3D, hitType, showerVertexPosition, showerSpineHitList))
            continue;
        */

        std::cout << "Finding shower start position..." << std::endl;
        // Find the 2D shower start position
        CartesianVector showerStartPosition(0.f, 0.f, 0.f);
        CartesianVector showerStartDirection(0.f, 0.f, 0.f);

        if (m_pShowerStartFinderTool->Run(pShowerPfo, peakDirection, hitType, showerSpineHitList, showerStartPosition, showerStartDirection) != STATUS_CODE_SUCCESS)
        {
            std::cout << "no shower start found" << std::endl;
            continue;
        }

        std::cout << "Finding shower start position..." << std::endl;
        CartesianVector cheatedShowerStartPosition(0.f, 0.f, 0.f);
        if (!this->GetShowerStart(pShowerPfo, hitType, showerSpineHitList, cheatedShowerStartPosition))
        {
            std::cout << "no cheated shower start found" << std::endl;
            continue;
        }

        ////////////////////////////////////////////
        PfoList visualisationPfo;
        visualisationPfo.push_back(pShowerPfo);
        PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &visualisationPfo, "ShowerPfos", RED);
        PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &showerSpineHitList, "ShowerSpine", RED);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &showerStartPosition, "showerStartPosition", BLACK, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &cheatedShowerStartPosition, "cheatedShowerStartPosition", BLACK, 2);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        break;
        ////////////////////////////////////////////
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MLShowerInitialRegionRefinementAlgorithm::GetShowerStart(const ParticleFlowObject *const pPfo, const HitType hitType, 
    const CaloHitList &showerSpineHitList, CartesianVector &showerStart) const
{
    // Combine HitLists
    CaloHitList combinedHitList;
    LArPfoHelper::GetCaloHits(pPfo, hitType, combinedHitList);

    for (const CaloHit *const pCaloHit : showerSpineHitList)
    {
        if (std::find(combinedHitList.begin(), combinedHitList.end(), pCaloHit) == combinedHitList.end())
            combinedHitList.push_back(pCaloHit);
    }

    // Get MCParticle List
    const MCParticleList *pMCParticleList(nullptr);

    if (PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList) != STATUS_CODE_SUCCESS)
        return false;

    if (!pMCParticleList || pMCParticleList->empty())
    {
        std::cout << "MLShowerInitialRegionRefinementAlgorithm: unable to find mc particle list " << m_mcParticleListName << std::endl;
        return false;
    }

    return (m_pCheatingShowerStartFinderTool->Run(pPfo, pMCParticleList, combinedHitList, hitType, showerStart) == STATUS_CODE_SUCCESS);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLShowerInitialRegionRefinementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoVertexListName", m_neutrinoVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameU", m_caloHitListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameV", m_caloHitListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameW", m_caloHitListNameW));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinShowerHits3D", m_minShowerHits3D));

    AlgorithmTool *pAlgorithmTool1(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "ShowerPeakDirectionFinder", pAlgorithmTool1));
    m_pShowerPeakDirectionFinderTool = dynamic_cast<PeakDirectionFinderTool *>(pAlgorithmTool1);

    if (!m_pShowerPeakDirectionFinderTool)
        return STATUS_CODE_INVALID_PARAMETER;

    AlgorithmTool *pAlgorithmTool2(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "ShowerSpineFinder", pAlgorithmTool2));
    m_pShowerSpineFinderTool = dynamic_cast<ShowerSpineFinderTool *>(pAlgorithmTool2);

    if (!m_pShowerSpineFinderTool)
        return STATUS_CODE_INVALID_PARAMETER;

    AlgorithmTool *pAlgorithmTool3(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "CheatingShowerStartFinder", pAlgorithmTool3));
    m_pCheatingShowerStartFinderTool = dynamic_cast<CheatingShowerStartFinderTool *>(pAlgorithmTool3);

    if (!m_pCheatingShowerStartFinderTool)
        return STATUS_CODE_INVALID_PARAMETER;

    AlgorithmTool *pAlgorithmTool4(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "ShowerStartFinder", pAlgorithmTool4));
    m_pShowerStartFinderTool = dynamic_cast<ShowerStartFinderTool *>(pAlgorithmTool4);

    if (!m_pShowerSpineFinderTool)
        return STATUS_CODE_INVALID_PARAMETER;

    return STATUS_CODE_SUCCESS;
}




//------------------------------------------------------------------------------------------------------------------------------------------
// List Manipulation
//------------------------------------------------------------------------------------------------------------------------------------------

void MLShowerInitialRegionRefinementAlgorithm::FillShowerPfoVector(PfoVector &showerPfoVector) const
{
    const PfoList *pPfoList(nullptr);

    if (PandoraContentApi::GetList(*this, m_showerPfoListName, pPfoList) != STATUS_CODE_SUCCESS)
        return;

    if (!pPfoList || pPfoList->empty())
    {
        std::cout << "MLShowerInitialRegionRefinementAlgorithm: unable to find shower pfo list " << m_showerPfoListName << std::endl;
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

StatusCode MLShowerInitialRegionRefinementAlgorithm::GetNeutrinoVertex(CartesianVector &nuVertex3D) const
{
    const VertexList *pNuVertexList(nullptr);
    const StatusCode statusCode(PandoraContentApi::GetList(*this, m_neutrinoVertexListName, pNuVertexList));

    if (statusCode != STATUS_CODE_SUCCESS)
        return statusCode;

    if (!pNuVertexList || (pNuVertexList->size() != 1))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "MLShowerInitialRegionRefinementAlgorithm: unable to find vertex list " << m_neutrinoVertexListName
                      << " if it does exist, it may have more than one nu vertex" << std::endl;

        return STATUS_CODE_NOT_INITIALIZED;
    }

    nuVertex3D = pNuVertexList->front()->GetPosition();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MLShowerInitialRegionRefinementAlgorithm::GetHitListOfType(const HitType hitType, const CaloHitList *&pCaloHitList) const
{
    const std::string typeHitListName(hitType == TPC_VIEW_U ? m_caloHitListNameU : hitType == TPC_VIEW_V ? m_caloHitListNameV : m_caloHitListNameW);

    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, typeHitListName, pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ShowerStartRefinementBaseTool: unable to find calo hit list " << typeHitListName << std::endl;

        return STATUS_CODE_NOT_INITIALIZED;
    }

    return STATUS_CODE_SUCCESS;
}


} // namespace lar_content
