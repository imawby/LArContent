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
    m_selectRecoHits{false},
    m_maxTierToExamine(4)
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
    //matchInfo.Print(mcHierarchy);

    MCParticleList nuParticles;
    mcHierarchy.GetRootMCParticles(nuParticles);

    // Try the fill tree function
    const LArHierarchyHelper::MCMatchesVector mcMatches(matchInfo.GetMatches(nuParticles.front()));

    for (const LArHierarchyHelper::MCMatches &mcMatch : mcMatches)
    {
        const LArHierarchyHelper::MCHierarchy::Node *const pPrimaryNode(mcMatch.GetMC());
        
        if (pPrimaryNode->GetHierarchyTier() != 1)
            continue;
        
        // Should it have been reco'd, and has it been?
        if ((!pPrimaryNode->IsReconstructable()) || (mcMatch.GetRecoMatches().empty()))
            continue;

        // Fill tree        
        LArHierarchyHelper::MCHierarchy::NodeVector hierarchyNodes;
        this->FoldMCNodes(pPrimaryNode, hierarchyNodes);  
        
        this->FillTree(pPrimaryNode, pPrimaryNode->GetMCParticles().front(), mcMatch.GetRecoMatches().front(), hierarchyNodes, mcMatches, 3);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SecondaryValidationAlgorithm::FoldMCNodes(const LArHierarchyHelper::MCHierarchy::Node *const pRootMCNode, LArHierarchyHelper::MCHierarchy::NodeVector &nodeVector)  
{
    nodeVector.push_back(pRootMCNode);

    for (const LArHierarchyHelper::MCHierarchy::Node *const pChildNode : pRootMCNode->GetChildren())
        this->FoldMCNodes(pChildNode, nodeVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SecondaryValidationAlgorithm::FillTree(const LArHierarchyHelper::MCHierarchy::Node *const pParentMCNode,
    const MCParticle *const pMCParent, const LArHierarchyHelper::RecoHierarchy::Node *const pParentRecoNode,
    const LArHierarchyHelper::MCHierarchy::NodeVector &hierarchyNodes,
    const LArHierarchyHelper::MCMatchesVector mcMatches, const int tierToExamine)
{
    if (pMCParent->GetDaughterList().empty() || (tierToExamine > m_maxTierToExamine))
        return;

    for (const MCParticle *const pMCChild : pMCParent->GetDaughterList())
    {        
        // Find MC node
        const LArHierarchyHelper::MCHierarchy::Node *pChildMCNode(nullptr);
        this->FindMCNode(pMCChild, hierarchyNodes, pChildMCNode);
        
        // If we can't find a node or if it isn't reco'able just skip it
        if ((pChildMCNode == nullptr) || (!pChildMCNode->IsReconstructable()))
        {   
            this->FillTree(pParentMCNode, pMCChild, pParentRecoNode, hierarchyNodes, mcMatches, tierToExamine);
            continue;
        }

        // If failed to reco an upstream parent
        if ((pParentMCNode == nullptr) || (pParentRecoNode == nullptr))
        {
            this->FillNullEntry(pMCParent, pMCChild, tierToExamine);
            this->FillTree(nullptr, pMCChild, nullptr, hierarchyNodes, mcMatches, (tierToExamine + 1));
            continue;
        }        
        else
        {        
            std::cout << "YES RECO: " << pChildMCNode->GetParticleId() << " GEN: " << tierToExamine << std::endl;

            const LArHierarchyHelper::MCMatches *pChildMatch(nullptr);
            this->FindMatch(pChildMCNode->GetMCParticles().front(), mcMatches, pChildMatch);
            this->FillEntry(pParentMCNode, pChildMCNode,  pChildMatch, mcMatches, tierToExamine);            
            this->FillTree(pChildMCNode, pMCChild, pChildMatch->GetRecoMatches().empty() ? nullptr : pChildMatch->GetRecoMatches().front(), hierarchyNodes, mcMatches, (tierToExamine + 1));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SecondaryValidationAlgorithm::FillNullEntry(const MCParticle *const pMCParent, const MCParticle *const pMCChild, const int hierarchyTier)
{
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ParentPDG", pMCParent->GetParticleId()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildPDG", pMCChild->GetParticleId()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "IsParentReco", 0));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildGeneration", hierarchyTier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildNHitsInHierarchyFrac", -999.f));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildCompleteness", -1.f));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildPurity", -1.f));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildNMatches", -999));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildTrueHitsU", -999));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildTrueHitsV", -999));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildTrueHitsW", -999));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildRecoHitsU", -999));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildRecoHitsV", -999));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildRecoHitsW", -999));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void SecondaryValidationAlgorithm::FillEntry(const LArHierarchyHelper::MCHierarchy::Node *const pParentMCNode,
    const LArHierarchyHelper::MCHierarchy::Node *const pChildMCNode, const LArHierarchyHelper::MCMatches *const pChildMatch, 
    const LArHierarchyHelper::MCMatchesVector &matchesVector, const int hierarchyTier)
{
    // Define variables to fill
    const int childNMatches(pChildMatch->GetRecoMatches().size());
    float childCompleteness(-1.f), childPurity(-1.f);
    int childTrueHitsU(-1), childTrueHitsV(-1), childTrueHitsW(-1);
    int childRecoHitsU(-1), childRecoHitsV(-1), childRecoHitsW(-1);

    // Fill reco node vars
    if (childNMatches != 0)
    {
        const LArHierarchyHelper::RecoHierarchy::Node *const pChildRecoNode(pChildMatch->GetRecoMatches().front());

        childCompleteness = pChildMatch->GetCompleteness(pChildRecoNode);
        childPurity = pChildMatch->GetPurity(pChildRecoNode);
        childTrueHitsU = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pChildMCNode->GetCaloHits());
        childTrueHitsV = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pChildMCNode->GetCaloHits());
        childTrueHitsW = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pChildMCNode->GetCaloHits());
        childRecoHitsU = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pChildRecoNode->GetCaloHits());
        childRecoHitsV = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pChildRecoNode->GetCaloHits());
        childRecoHitsW = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pChildRecoNode->GetCaloHits());
    }

    const int nHitsInHierarchy(this->GetHitsInUpstreamHierarchy(pChildMCNode, matchesVector));
    const float hitsInHierarchyFrac(static_cast<float>(nHitsInHierarchy) / pChildMCNode->GetCaloHits().size());

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ParentPDG", pParentMCNode->GetParticleId()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildPDG", pChildMCNode->GetParticleId()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "IsParentReco", 1));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildGeneration", hierarchyTier));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildNHitsInHierarchyFrac", hitsInHierarchyFrac));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildCompleteness", childCompleteness));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildPurity", childPurity));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildNMatches", childNMatches));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildTrueHitsU", childTrueHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildTrueHitsV", childTrueHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildTrueHitsW", childTrueHitsW));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildRecoHitsU", childRecoHitsU));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildRecoHitsV", childRecoHitsV));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ChildRecoHitsW", childRecoHitsW));
    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
}

//------------------------------------------------------------------------------------------------------------------------------------------

int SecondaryValidationAlgorithm::GetHitsInUpstreamHierarchy(const LArHierarchyHelper::MCHierarchy::Node *const pRootMCNode, 
    const LArHierarchyHelper::MCMatchesVector &matchesVector)
{
    int hitsInHierarchy(0); 

    MCParticleList ancestorMCParticles;
    LArMCParticleHelper::GetAllAncestorMCParticles(pRootMCNode->GetMCParticles().front(), ancestorMCParticles);

    for (const MCParticle *const pMCAncestor : ancestorMCParticles)
    {
        const LArHierarchyHelper::MCMatches *pRecoMatch(nullptr);

        if (!this->FindMatch(pMCAncestor, matchesVector, pRecoMatch))
            continue;

        for (const LArHierarchyHelper::RecoHierarchy::Node *const pRecoNode : pRecoMatch->GetRecoMatches())
        {
            CaloHitVector intersection;
            std::set_intersection(pRecoNode->GetCaloHits().begin(), pRecoNode->GetCaloHits().end(), 
                                  pRootMCNode->GetCaloHits().begin(), pRootMCNode->GetCaloHits().end(), std::back_inserter(intersection));

            hitsInHierarchy += intersection.size();
        }
    }

    return hitsInHierarchy;
}

//------------------------------------------------------------------------------------------------------------------------------------------
// Find functions
//------------------------------------------------------------------------------------------------------------------------------------------

bool SecondaryValidationAlgorithm::FindMCNode(const MCParticle *const pMCParticle, const LArHierarchyHelper::MCHierarchy::NodeVector &hierarchyMCNodes, 
    const LArHierarchyHelper::MCHierarchy::Node *&pMCNode)
{
    for (const LArHierarchyHelper::MCHierarchy::Node *const pNode : hierarchyMCNodes)
    {
        if (pNode->GetMCParticles().front() == pMCParticle)
        {
            pMCNode = pNode;
            return true;
        }
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool SecondaryValidationAlgorithm::FindMatch(const MCParticle *const pMCParticle, const LArHierarchyHelper::MCMatchesVector &matchesVector, 
    const LArHierarchyHelper::MCMatches *&pRecoMatch)  
{
    // Find reco node
    for (const LArHierarchyHelper::MCMatches &mcMatch : matchesVector)
    {
        if (mcMatch.GetMC()->GetMCParticles().front() == pMCParticle)
        {
            pRecoMatch = &mcMatch;
            return true;
        }
    }

    return false;
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
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxTierToExamine", m_maxTierToExamine));    

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
