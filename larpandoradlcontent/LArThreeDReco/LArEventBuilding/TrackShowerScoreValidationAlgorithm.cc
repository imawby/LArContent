/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/TrackShowerScoreValidationAlgorithm.cc
 *
 *  @brief  Implementation of the MLP neutrino hierarchy algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoradlcontent/LArCheating/MLPCheatHierarchyTool.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/TrackShowerScoreValidationAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{


//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

TrackShowerScoreValidationAlgorithm::TrackShowerScoreValidationAlgorithm() :
    m_eventID(0),
    m_pfoListNames({"TrackParticles3D", "ShowerParticles3D"}),
    m_validationFileName("TrackShowerScoreValidationFile.root"),
    m_validationTreeName("TrackShowerTree")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackShowerScoreValidationAlgorithm::~TrackShowerScoreValidationAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_validationTreeName.c_str(), m_validationFileName.c_str(), "UPDATE"));
    }
    catch (...) {}
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackShowerScoreValidationAlgorithm::Run()
{
    // Do ID stuff
    ++m_eventID;
        
    // Get truth
    PfoToMCParticleMap matchingMap; ChildToParentPfoMap childToParentPfoMap;
    m_cheatHierarchyTool->FillHierarchyMap(this, matchingMap, childToParentPfoMap);

    // Get All pfos
    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);

        if (PandoraContentApi::GetList(*this, pfoListName, pPfoList) != STATUS_CODE_SUCCESS)
            continue;

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {

            if (matchingMap.find(pPfo) == matchingMap.end())
                continue;
            
            const MCParticle *const pMCParticle(matchingMap.at(pPfo));
            const int matchedPDG(pMCParticle->GetParticleId());

            float trackShowerScore(-999.f);
            const PropertiesMap &pfoMeta(pPfo->GetPropertiesMap());

            if (pfoMeta.find("TrackScore") != pfoMeta.end())
                trackShowerScore = pfoMeta.at("TrackScore");

            const int nSpacepoints(this->GetNSpacepoints(pPfo));
            
            // Fill trees
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_validationTreeName.c_str(), "EventID", m_eventID));            
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_validationTreeName.c_str(), "MatchedPDG", matchedPDG));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_validationTreeName.c_str(), "TrackShowerScore", trackShowerScore));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_validationTreeName.c_str(), "NSpacepoints", nSpacepoints));                        
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_validationTreeName.c_str()));
            
        }
    }
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float TrackShowerScoreValidationAlgorithm::GetNSpacepoints(const ParticleFlowObject *const pPfo) const
{
    ClusterList clusterList3D;
    LArPfoHelper::GetThreeDClusterList(pPfo, clusterList3D);

    int total3DHits(0);

    for (const Cluster *const pCluster3D : clusterList3D)
        total3DHits += pCluster3D->GetNCaloHits();

    return total3DHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------
StatusCode TrackShowerScoreValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ValidationFileName", m_validationFileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ValidationTreeName", m_validationTreeName));
    
    AlgorithmTool *pAlgorithmTool(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "MLPCheatHierarchyTool", pAlgorithmTool));
    m_cheatHierarchyTool = dynamic_cast<MLPCheatHierarchyTool *>(pAlgorithmTool);

    if (!m_cheatHierarchyTool)
        return STATUS_CODE_INVALID_PARAMETER;

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
