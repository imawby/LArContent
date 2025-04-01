/**
 *  @file   larpandoracontent/LArMonitoring/SecondaryValidationAlgorithm.cc
 *
 *  @brief  Implementation of the secondary validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArHierarchyHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "larpandoracontent/LArMonitoring/SecondaryValidationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

SecondaryValidationAlgorithm::SecondaryValidationAlgorithm() :
    m_writeFile(true),
    m_fileName("HierarchyRecoPerformance.root"),
    m_treeName("tree"),
    m_minPurity{0.8f},
    m_minCompleteness{0.65f},
    m_minRecoHits{30},
    m_minRecoHitsPerView{10},
    m_minRecoGoodViews{2},
    m_removeRecoNeutrons{true},
    m_selectRecoHits{false}
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

SecondaryValidationAlgorithm::~SecondaryValidationAlgorithm()
{
    if (m_writeFile)
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName.c_str(), m_fileName.c_str(), "UPDATE"));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SecondaryValidationAlgorithm::Run()
{
    //++m_event;
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
    const PfoList *pPfoList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    LArHierarchyHelper::FoldingParameters foldParameters;
    foldParameters.m_foldToLeadingShowers = true;

    const LArHierarchyHelper::MCHierarchy::ReconstructabilityCriteria recoCriteria(
        m_minRecoHits, m_minRecoHitsPerView, m_minRecoGoodViews, m_removeRecoNeutrons);

    LArHierarchyHelper::MCHierarchy mcHierarchy(recoCriteria);
    LArHierarchyHelper::FillMCHierarchy(*pMCParticleList, *pCaloHitList, foldParameters, mcHierarchy);
    LArHierarchyHelper::RecoHierarchy recoHierarchy;
    LArHierarchyHelper::FillRecoHierarchy(*pPfoList, foldParameters, recoHierarchy);
    const LArHierarchyHelper::QualityCuts quality(m_minPurity, m_minCompleteness, m_selectRecoHits);
    LArHierarchyHelper::MatchInfo matchInfo(mcHierarchy, recoHierarchy, quality);
    LArHierarchyHelper::MatchHierarchies(matchInfo);
    matchInfo.Print(mcHierarchy);

    MCParticleList nuParticles;
    mcHierarchy.GetRootMCParticles(nuParticles);

    // Loop through primaries
    if (!nuParticles.empty())
    {
        const LArHierarchyHelper::MCHierarchy::NodeVector primaryNodes(mcHierarchy.GetInteractions(nuParticles.front()));
        const LArHierarchyHelper::MCMatchesVector mcMatches(matchInfo.GetMatches(nuParticles.front()));

        for (const LArHierarchyHelper::MCMatches &mcMatch : mcMatches)
        {
            // Has it been reco'd?
            if (mcMatch.GetRecoMatches().empty())
                continue;

            // Define vectors to save into...
            std::vector<int> parentPDG;
            std::vector<int> childPDG;
            std::vector<int> childGen;
            std::vector<float> childHitsInHierarchyFrac;
            std::vector<int> childNMatches;
            std::vector<float> childCompleteness;
            std::vector<float> childPurity;
            std::vector<int> childTrueHitsU, childTrueHitsV, childTrueHitsW;
            std::vector<int> childRecoHitsU, childRecoHitsV, childRecoHitsW;

            const LArHierarchyHelper::MCHierarchy::Node *const pPrimaryNode(mcMatch.GetMC());
            const LArHierarchyHelper::MCHierarchy::NodeVector &childNodes(pPrimaryNode->GetChildren());

            std::cout << "PDG: " << pPrimaryNode->GetParticleId() << std::endl;

            for (const LArHierarchyHelper::MCHierarchy::Node *const pChildTrueNode : childNodes)
            {
                if (!pChildTrueNode->IsReconstructable())
                    continue;

                int thisParentPDG(pPrimaryNode->GetParticleId());
                int thisChildPDG(pChildTrueNode->GetParticleId());
                int thisChildGen(3);
                float thisChildHitsInHierarchyFrac(-1.f);
                int thisChildNMatches(-1);
                float thisChildCompleteness(-1.f);
                float thisChildPurity(-1.f);
                int thisChildTrueHitsU(-1), thisChildTrueHitsV(-1), thisChildTrueHitsW(-1);
                int thisChildRecoHitsU(-1), thisChildRecoHitsV(-1), thisChildRecoHitsW(-1);

                // Find child reco node (if it has one!)
                for (const LArHierarchyHelper::MCMatches &mcMatchTemp : mcMatches)
                {
                    if (mcMatchTemp.GetMC() == pChildTrueNode)
                    {
                        thisChildNMatches = mcMatchTemp.GetRecoMatches().size();

                        if (thisChildNMatches != 0)
                        {
                            // Get the best match
                            const LArHierarchyHelper::RecoHierarchy::Node *const pChildRecoNode(mcMatchTemp.GetRecoMatches().front());
                            thisChildCompleteness = mcMatchTemp.GetCompleteness(pChildRecoNode);
                            thisChildPurity = mcMatchTemp.GetPurity(pChildRecoNode);
                            thisChildTrueHitsU = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pChildTrueNode->GetCaloHits());
                            thisChildTrueHitsV = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pChildTrueNode->GetCaloHits());
                            thisChildTrueHitsW = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pChildTrueNode->GetCaloHits());
                            thisChildRecoHitsU = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pChildRecoNode->GetCaloHits());
                            thisChildRecoHitsV = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pChildRecoNode->GetCaloHits());
                            thisChildRecoHitsW = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pChildRecoNode->GetCaloHits());
                        }

                        break;
                    }
                }

                // Get hits in hierarchy?
                int thisChildHitsInHierarchy(0);
                for (const LArHierarchyHelper::RecoHierarchy::Node *const pParentRecoNode : mcMatch.GetRecoMatches())
                {
                    CaloHitVector intersection;
                    std::set_intersection(pParentRecoNode->GetCaloHits().begin(), pParentRecoNode->GetCaloHits().end(), 
                                          pChildTrueNode->GetCaloHits().begin(), pChildTrueNode->GetCaloHits().end(), std::back_inserter(intersection));

                    thisChildHitsInHierarchy += intersection.size();
                }

                if (!pChildTrueNode->GetCaloHits().empty())
                    thisChildHitsInHierarchyFrac = static_cast<float>(thisChildHitsInHierarchy) / pChildTrueNode->GetCaloHits().size(); 

                std::cout << " ----- PDG: " << pChildTrueNode->GetParticleId() << std::endl;

                parentPDG.push_back(thisParentPDG);
                childPDG.push_back(thisChildPDG);
                childGen.push_back(thisChildGen);
                childHitsInHierarchyFrac.push_back(thisChildHitsInHierarchyFrac);
                childNMatches.push_back(thisChildNMatches);
                childCompleteness.push_back(thisChildCompleteness);
                childPurity.push_back(thisChildPurity);
                childTrueHitsU.push_back(thisChildTrueHitsU);
                childTrueHitsV.push_back(thisChildTrueHitsV);
                childTrueHitsW.push_back(thisChildTrueHitsW);
                childRecoHitsU.push_back(thisChildRecoHitsU);
                childRecoHitsV.push_back(thisChildRecoHitsV);
                childRecoHitsW.push_back(thisChildRecoHitsW);
            }

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "PrimaryPDG", &parentPDG));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecondaryPDG", &childPDG));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecondaryNHitsInParent", &childHitsInHierarchyFrac));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecondaryCompleteness", &childCompleteness));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecondaryPurity", &childPurity));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecondaryNMatches", &childNMatches));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecondaryTrueHitsU", &childTrueHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecondaryTrueHitsV", &childTrueHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecondaryTrueHitsW", &childTrueHitsW));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecondaryRecoHitsU", &childRecoHitsU));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecondaryRecoHitsV", &childRecoHitsV));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "SecondaryRecoHitsW", &childRecoHitsW));
            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode SecondaryValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    if (m_caloHitListName.empty())
        m_caloHitListName = "CaloHitList2D";
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));
    if (m_pfoListName.empty())
        m_pfoListName = "RecreatedPfos";

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteFile", m_writeFile));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinPurity", m_minPurity));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCompleteness", m_minCompleteness));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoHits", m_minRecoHits));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoHitsPerView", m_minRecoHitsPerView));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinRecoGoodViews", m_minRecoGoodViews));
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RemoveRecoNeutrons", m_removeRecoNeutrons));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
