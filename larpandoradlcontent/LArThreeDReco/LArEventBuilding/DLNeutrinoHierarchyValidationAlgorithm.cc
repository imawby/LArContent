/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLNeutrinoHierarchyValidationAlgorithm.cc
 *
 *  @brief  Implementation of the DL neutrino hierarchy algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoradlcontent/LArCheating/DLCheatHierarchyTool.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"

#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLNeutrinoHierarchyValidationAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DLNeutrinoHierarchyValidationAlgorithm::DLNeutrinoHierarchyValidationAlgorithm() :
    m_validationFileName("DLHierarchyValidationFile.root"),
    m_validationTreeName("ValidationTree"),    
    m_mcParticleListName("Input"),
    m_neutrinoPfoListName("NeutrinoParticles3D"),
    m_pfoListNames({"TrackParticles3D", "ShowerParticles3D"}),
    m_minClusterSize(5),
    m_slidingFitWindow(20)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DLNeutrinoHierarchyValidationAlgorithm::~DLNeutrinoHierarchyValidationAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_validationTreeName.c_str(), m_validationFileName.c_str(), "UPDATE"));
    }
    catch (...) {}    
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLNeutrinoHierarchyValidationAlgorithm::Run()
{
    const ParticleFlowObject *pNeutrinoPfo(nullptr);

    if (!this->GetNeutrinoPfo(pNeutrinoPfo))
        return STATUS_CODE_SUCCESS;

    HierarchyPfoVector trackPfos, showerPfos;

    this->FillTrackShowerVectors(pNeutrinoPfo, trackPfos, showerPfos);
    this->Validate(pNeutrinoPfo, trackPfos, showerPfos);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DLNeutrinoHierarchyValidationAlgorithm::GetNeutrinoPfo(const ParticleFlowObject *&pNeutrinoPfo) const
{
    const PfoList *pPfoList(nullptr);

    PANDORA_THROW_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, m_neutrinoPfoListName, pPfoList));

    if (!pPfoList || pPfoList->empty())
        return false;

    // ATTN Enforces that only one pfo, of neutrino-type, be in the specified input list
    pNeutrinoPfo = (1 == pPfoList->size()) ? *(pPfoList->begin()) : nullptr;

    if (!pNeutrinoPfo || !LArPfoHelper::IsNeutrino(pNeutrinoPfo))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLNeutrinoHierarchyValidationAlgorithm::FillTrackShowerVectors(
    const ParticleFlowObject *const pNeutrinoPfo, HierarchyPfoVector &trackPfos, HierarchyPfoVector &showerPfos) const
{
    // Sliding fit shenanigans
    const float pitchU{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U)};
    const float pitchV{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V)};
    const float pitchW{LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W)};
    const float slidingFitPitch(std::max({pitchU, pitchV, pitchW}));

    // Get All pfos
    PfoVector pfoVector;
    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        if (PandoraContentApi::GetList(*this, pfoListName, pPfoList) != STATUS_CODE_SUCCESS)
            continue;

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            // Apply hit cut
            if (LArPfoHelper::GetNumberOfThreeDHits(pPfo) >= m_minClusterSize)
                pfoVector.emplace_back(pPfo);
        }
    }

    // Sort to maintain reproducibility
    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);

    for (const ParticleFlowObject *const pPfo : pfoVector)
    {
        // Attempt sliding linear fit
        try
        {
            ClusterList clusters3D;
            LArPfoHelper::GetThreeDClusterList(pPfo, clusters3D);

            if (clusters3D.size() != 1)
                continue;

            const ThreeDSlidingFitResult slidingFitResult(*clusters3D.begin(), m_slidingFitWindow, slidingFitPitch);

            // We need directions...
            ExtremalPoint upstreamPoint, downstreamPoint;

            // Create track/shower objects
            if (pPfo->GetParticleId() == 13)
            {
                trackPfos.push_back(HierarchyPfo(pPfo, slidingFitResult, upstreamPoint, downstreamPoint));
            }
            else if (pPfo->GetParticleId() == 11)
            {
                showerPfos.push_back(HierarchyPfo(pPfo, slidingFitResult, upstreamPoint, downstreamPoint));
            }
        }
        catch (...)
        {
            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLNeutrinoHierarchyValidationAlgorithm::Validate(const ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoVector &trackPfos,
    const HierarchyPfoVector &showerPfos) const
{
    // Get truth
    PfoToMCParticleMap matchingMap; ChildToParentPfoMap childToParentPfoMap;
    m_cheatHierarchyTool->FillHierarchyMap(this, matchingMap, childToParentPfoMap);

    const float nuVertexAccuracy(this->GetNuVertexAccuracy(pNeutrinoPfo));
    
    // Get All pfos
    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        if (PandoraContentApi::GetList(*this, pfoListName, pPfoList) != STATUS_CODE_SUCCESS)
            continue;

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            int isTrack(pPfo->GetParticleId() == 13);
            
            // Do we know the MCParticle?
            if (matchingMap.find(pPfo) == matchingMap.end())
                continue;
            
            // Is in matching list?
            if (childToParentPfoMap.find(pPfo) == childToParentPfoMap.end())
                continue;

            // Is one of our target PFPs?
            if ((std::find(trackPfos.begin(), trackPfos.end(), pPfo) == trackPfos.end()) && 
                (std::find(showerPfos.begin(), showerPfos.end(), pPfo) == showerPfos.end()))
            {
                continue;
            }
            
            // Get true info
            const MCParticle *const pMCParticle(matchingMap.at(pPfo));
            const int matchedPDG(pMCParticle->GetParticleId());
            const ParticleFlowObject *const pTrueParentPfo(childToParentPfoMap.at(pPfo).first);
            const int trueVisibleGen(childToParentPfoMap.at(pPfo).second);            
            
            PfoList parentList(pPfo->GetParentPfoList());

            if (parentList.size() != 1)
                continue;

            const ParticleFlowObject *const pRecoParent(*parentList.begin());
            
            // Is parent one of our target PFPs?
            if ((pTrueParentPfo != pNeutrinoPfo) && 
                (std::find(trackPfos.begin(), trackPfos.end(), pTrueParentPfo) == trackPfos.end()) && 
                (std::find(showerPfos.begin(), showerPfos.end(), pTrueParentPfo) == showerPfos.end()))
            {
                continue;
            }

            const int recoGen(LArPfoHelper::GetHierarchyTier(pPfo) + 1);
            const int hasCorrectParent(pRecoParent == pTrueParentPfo ? 1 : 0);

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_validationTreeName.c_str(), "IsTrack", isTrack));            
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_validationTreeName.c_str(), "TruePDG", matchedPDG));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_validationTreeName.c_str(), "TrueVisibleGen", trueVisibleGen));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_validationTreeName.c_str(), "HasCorrectParent", hasCorrectParent));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_validationTreeName.c_str(), "RecoGen", recoGen));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_validationTreeName.c_str(), "NuVertexAccuracy", nuVertexAccuracy));                        
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_validationTreeName.c_str()));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float DLNeutrinoHierarchyValidationAlgorithm::GetNuVertexAccuracy(const ParticleFlowObject *const pNeutrinoPfo) const
{
    const MCParticleList *pMCParticleList(nullptr);
    if (PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList) != STATUS_CODE_SUCCESS)
        return std::numeric_limits<float>::max();

    if (!pMCParticleList || pMCParticleList->empty())
        return std::numeric_limits<float>::max();        

    // Apply nu vertex accuracy requirement
    MCParticleVector mcNeutrinoVector;
    LArMCParticleHelper::GetTrueNeutrinos(pMCParticleList, mcNeutrinoVector);

    if (mcNeutrinoVector.size() != 1)
        return std::numeric_limits<float>::max();                

    const CartesianVector &trueNuVertex(mcNeutrinoVector.front()->GetVertex());
    const Vertex *const pNeutrinoVertex(LArPfoHelper::GetVertex(pNeutrinoPfo));
    const CartesianVector recoNuVertex(pNeutrinoVertex->GetPosition().GetX(), pNeutrinoVertex->GetPosition().GetY(), pNeutrinoVertex->GetPosition().GetZ());    

    return (trueNuVertex - recoNuVertex).GetMagnitudeSquared();
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLNeutrinoHierarchyValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
                                    STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterSize", m_minClusterSize));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoPfoListName", m_neutrinoPfoListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ValidationFileName", m_validationFileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ValidationTreeName", m_validationTreeName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterSize", m_minClusterSize));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    AlgorithmTool *pAlgorithmTool(nullptr);
    pAlgorithmTool = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "DLCheatHierarchyTool", pAlgorithmTool));
    m_cheatHierarchyTool = dynamic_cast<DLCheatHierarchyTool *>(pAlgorithmTool);

    if (!m_cheatHierarchyTool)
        return STATUS_CODE_INVALID_PARAMETER;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
