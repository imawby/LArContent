/**
 *  @file   larpandoracontent/LArHelpers/LArMuonLeadingHelper.cc
 *
 *  @brief  Implementation of the lar delta ray helper class.
 *
 *  $Log: $
 */

#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMuonLeadingHelper.h"

#include "Objects/ParticleFlowObject.h"

namespace lar_content
{

using namespace pandora;

LArMuonLeadingHelper::ValidationParameters::ValidationParameters() :
    LArMCParticleHelper::PrimaryParameters()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------    

bool LArMuonLeadingHelper::IsDeltaRay(const MCParticle *const pMCParticle)
{
    const MCParticleList parentList(pMCParticle->GetParentList());

    if(parentList.empty())
        return false;

    if (1 != parentList.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if (!LArMCParticleHelper::IsCosmicRay(parentList.front()))
        return false;

    const LArMCParticle *const pLArMCParticle(dynamic_cast<const LArMCParticle*>(pMCParticle));
    
    return pLArMCParticle->GetIsDR();
}

//------------------------------------------------------------------------------------------------------------------------------------------    

bool LArMuonLeadingHelper::IsMichel(const MCParticle *const pMCParticle)
{
    if (std::abs(pMCParticle->GetParticleId()) != 11)
        return false;
    
    const MCParticleList parentList(pMCParticle->GetParentList());

    if(parentList.empty())
        return false;
    
    if (1 != parentList.size())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    
    if (!LArMCParticleHelper::IsCosmicRay(parentList.front()))
        return false;    
    
    /*
    if (std::abs(parentList.front()->GetParticleId()) != 13)
        return false;
    */
    const LArMCParticle *const pLArMCParticle(dynamic_cast<const LArMCParticle*>(pMCParticle));

    return pLArMCParticle->GetIsDecay();
}

//------------------------------------------------------------------------------------------------------------------------------------------    

bool LArMuonLeadingHelper::IsLeading(const MCParticle *const pMCParticle)
{
    const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));
    
    return ((LArMCParticleHelper::GetHierarchyTier(pMCParticle) == 1) && (LArMCParticleHelper::IsCosmicRay(pParentMCParticle)));
    //return ((LArMCParticleHelper::GetHierarchyTier(pMCParticle) == 1) && (std::abs(pParentMCParticle->GetParticleId()) == 13));
}     

//------------------------------------------------------------------------------------------------------------------------------------------

const MCParticle *LArMuonLeadingHelper::GetLeadingParticle(const MCParticle *const pMCParticle)
{
    // Navigate upward through MC daughter/parent links - collect this particle and all its parents
    MCParticleVector mcParticleVector;

    const MCParticle *pParentMCParticle = pMCParticle;
    mcParticleVector.push_back(pParentMCParticle);

    while (!pParentMCParticle->GetParentList().empty())
    {
        if (1 != pParentMCParticle->GetParentList().size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pParentMCParticle = *(pParentMCParticle->GetParentList().begin());
        mcParticleVector.push_back(pParentMCParticle);
    }

    return *(mcParticleVector.end() - 2);
}    

//------------------------------------------------------------------------------------------------------------------------------------------    

void LArMuonLeadingHelper::GetMCToLeadingMap(const MCParticleList *const pMCParticleList, LArMCParticleHelper::MCRelationMap &mcToLeadingMap)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));

        if (pMCParticle == pParentMCParticle)
            continue;
        
        if (!LArMCParticleHelper::IsCosmicRay(pParentMCParticle))
            continue;
        
        /*
        if (std::abs(pParentMCParticle->GetParticleId()) != 13)
            continue;
        */      
        const MCParticle *const pLeadingMCParticle(LArMuonLeadingHelper::GetLeadingParticle(pMCParticle));
            
        mcToLeadingMap[pMCParticle] = pLeadingMCParticle;
    }
}    

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::RemoveMuonPfosFromList(const PfoList *const pPfoList, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &unfoldedInterepretedMatchingMap,
    PfoList &outputList)
{
    MCParticleVector reconstructableMCParticles;
    for (auto &entry : unfoldedInterepretedMatchingMap)
        reconstructableMCParticles.push_back(entry.first);

    std::sort(reconstructableMCParticles.begin(), reconstructableMCParticles.end(), LArMCParticleHelper::SortByMomentum);

    PfoList muonPfos;
    for (const MCParticle *const pMCParticle : reconstructableMCParticles)
    {
        if (std::fabs(pMCParticle->GetParticleId()) != 13)
            continue;

        const auto iter(unfoldedInterepretedMatchingMap.find(pMCParticle));

        for (auto &matchedPfoHitsPair : iter->second)
            muonPfos.push_back(matchedPfoHitsPair.first);
    }

    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        if (std::find(muonPfos.begin(), muonPfos.end(), pPfo) == muonPfos.end())
            outputList.push_back(pPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::SelectNonMuonLeadingPfos(const PfoList &inputPfoList, PfoList &outputList)
{    
    for (const ParticleFlowObject *const pPfo : inputPfoList)
    {
        const PfoList &parentPfoList(pPfo->GetParentPfoList());

        if (parentPfoList.size() > 1)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        if (parentPfoList.empty())
        {
            outputList.push_back(pPfo); 
        }
        else
        {
            if (std::find(inputPfoList.begin(), inputPfoList.end(), parentPfoList.front()) == inputPfoList.end())
                outputList.push_back(pPfo);
        }
    }
}    

//------------------------------------------------------------------------------------------------------------------------------------------    

void LArMuonLeadingHelper::SelectReconstructableLeadingParticles(const MCParticleList *pMCParticleList, const CaloHitList *pCaloHitList, const ValidationParameters &parameters,
    LArMCParticleHelper::MCContributionMap &selectedMCParticlesToHitsMap)
{
    // Obtain map: [mc particle -> leading muon child mc]
    LArMCParticleHelper::MCRelationMap mcToLeadingMCMap;
    LArMuonLeadingHelper::GetMCToLeadingMap(pMCParticleList, mcToLeadingMCMap);

    // Select reconstructable hits, e.g. remove those downstream of a neutron
    // Unless selectInputHits == false
    CaloHitList selectedCaloHitList;
    LArMCParticleHelper::SelectCaloHits(pCaloHitList, mcToLeadingMCMap, selectedCaloHitList, parameters.m_selectInputHits, parameters.m_maxPhotonPropagation);//std::numeric_limits<float>::max());

    // Obtain maps: [hit -> leading muon child mc], [leading muon child mc -> list of hits]
    LArMCParticleHelper::CaloHitToMCMap trueHitToLeadingMCMap;
    LArMCParticleHelper::MCContributionMap leadingMCToTrueHitListMap;
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(&selectedCaloHitList, mcToLeadingMCMap, trueHitToLeadingMCMap, leadingMCToTrueHitListMap);

    // Obtain vector: all mc particles
    MCParticleVector leadingMCVector;
    LArMuonLeadingHelper::SelectLeadingMCParticles(pMCParticleList, leadingMCVector);
    
    // Ensure the MCParticles have enough "good" hits to be reconstructed
    LArMCParticleHelper::SelectParticlesByHitCount(leadingMCVector, leadingMCToTrueHitListMap, mcToLeadingMCMap, parameters, selectedMCParticlesToHitsMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::GetPfoMatchContamination(const MCParticle *const pLeadingParticle, const CaloHitList &matchedPfoHitList,
    CaloHitList &parentTrackHits, CaloHitList &otherTrackHits, CaloHitList &otherShowerHits)
{
    const MCParticle *const pParentCosmicRay(LArMCParticleHelper::GetParentMCParticle(pLeadingParticle));
    
    for (const CaloHit *const pCaloHit : matchedPfoHitList)
    {
        const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

        if (LArMCParticleHelper::IsCosmicRay(pHitParticle))
        {
            if (pHitParticle == pParentCosmicRay)
            {
                parentTrackHits.push_back(pCaloHit);
            }
            else
            {
                otherTrackHits.push_back(pCaloHit);
            }
        }
        else
        {
            const MCParticle *const pHitLeadingParticle(LArMuonLeadingHelper::GetLeadingParticle(pHitParticle));

            if (pHitLeadingParticle != pLeadingParticle)
                otherShowerHits.push_back(pCaloHit);
        }
    }
} 

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::GetMuonPfoContaminationContribution(const CaloHitList &cosmicRayPfoHitList, const CaloHitList &leadingMCHitList,
    CaloHitList &leadingHitsInParentCosmicRay)
{
    for (const CaloHit *const pCaloHit : cosmicRayPfoHitList)
    {
        if (std::find(leadingMCHitList.begin(), leadingMCHitList.end(), pCaloHit) != leadingMCHitList.end())
            leadingHitsInParentCosmicRay.push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::SelectLeadingMCParticles(const MCParticleList *pMCParticleList, MCParticleVector &selectedParticles)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (LArMuonLeadingHelper::IsLeading(pMCParticle))
            selectedParticles.push_back(pMCParticle);
    }

    std::sort(selectedParticles.begin(), selectedParticles.end(), LArMCParticleHelper::SortByMomentum);
}
    
} // namespace lar_content
