/**
 *  @file   larpandoracontent/LArHelpers/LArMuonLeadingHelper.cc
 *
 *  @brief  Implementation of the lar delta ray helper class.
 *
 *  $Log: $
 */

#include "Helpers/MCParticleHelper.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArMuonLeadingHelper.h"

#include "Pandora/PdgTable.h"

#include "PandoraMonitoringApi.h"

#include "Objects/ParticleFlowObject.h"
#include "Objects/CaloHit.h"

namespace lar_content
{

using namespace pandora;

LArMuonLeadingHelper::ValidationParameters::ValidationParameters() :
    LArMCParticleHelper::PrimaryParameters(),
    m_maxBremsstrahlungSeparation(std::numeric_limits<float>::max())
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
    
    const LArMCParticle *const pLArMCParticle(dynamic_cast<const LArMCParticle*>(pMCParticle));

    return pLArMCParticle->GetIsDecay();
}

//------------------------------------------------------------------------------------------------------------------------------------------    

bool LArMuonLeadingHelper::IsLeading(const MCParticle *const pMCParticle)
{
    const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));
    
    return ((LArMCParticleHelper::GetHierarchyTier(pMCParticle) == 1) && (LArMCParticleHelper::IsCosmicRay(pParentMCParticle)));
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

        if (!LArMCParticleHelper::IsCosmicRay(pParentMCParticle))
            continue;

        // For the CRs: fold hits to themselves, for the DRs: fold hits to the leading non-muon MC
        if (pMCParticle == pParentMCParticle)
        {
            mcToLeadingMap[pMCParticle] = pMCParticle;
        }
        else
        {
            const MCParticle *const pLeadingMCParticle(LArMuonLeadingHelper::GetLeadingParticle(pMCParticle));
            mcToLeadingMap[pMCParticle] = pLeadingMCParticle;
        }
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
    LArMCParticleHelper::MCContributionMap &selectedMCParticlesToHitsMap, const Pandora &pandora)
{
    // Obtain map: [mc particle -> leading muon child mc]
    LArMCParticleHelper::MCRelationMap mcToLeadingMCMap;
    LArMuonLeadingHelper::GetMCToLeadingMap(pMCParticleList, mcToLeadingMCMap);

    // Select reconstructable hits, e.g. remove those downstream of a neutron
    // Unless selectInputHits == false
    CaloHitList selectedCaloHitList;
    LeadingMCParticleToPostPhotonHitLists leadingMCParticleToPostPhotonHitLists;
    LArMuonLeadingHelper::SelectCaloHits(pCaloHitList, mcToLeadingMCMap, selectedCaloHitList, parameters.m_selectInputHits, parameters.m_minHitSharingFraction, leadingMCParticleToPostPhotonHitLists);

    // Obtain maps: [hit -> leading muon child mc], [leading muon child mc -> list of hits]
    LArMCParticleHelper::CaloHitToMCMap trueHitToLeadingMCMap;
    LArMCParticleHelper::MCContributionMap leadingMCToTrueHitListMap;
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(&selectedCaloHitList, mcToLeadingMCMap, trueHitToLeadingMCMap, leadingMCToTrueHitListMap);

    // Add back in post bremsstrahlung hits that are close to the non-muon leading
    LArMuonLeadingHelper::AddInReconstructablePostPhotonHits(leadingMCParticleToPostPhotonHitLists, parameters.m_maxBremsstrahlungSeparation, leadingMCToTrueHitListMap, pandora);

    // Obtain vector: all mc particles
    MCParticleVector leadingMCVector;
    LArMuonLeadingHelper::SelectLeadingMCParticles(pMCParticleList, leadingMCVector);
    
    // Ensure the MCParticles have enough "good" hits to be reconstructed
    LArMCParticleHelper::SelectParticlesByHitCount(leadingMCVector, leadingMCToTrueHitListMap, mcToLeadingMCMap, parameters, selectedMCParticlesToHitsMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::SelectCaloHits(const CaloHitList *const pCaloHitList, const LArMCParticleHelper::MCRelationMap &mcToTargetMCMap,
    CaloHitList &selectedCaloHitList, const bool selectInputHits, const float minHitSharingFraction, LeadingMCParticleToPostPhotonHitLists &leadingMCParticleToPostPhotonHitLists)
{
    if (!selectInputHits)
    {
        selectedCaloHitList.insert(selectedCaloHitList.end(), pCaloHitList->begin(), pCaloHitList->end());
        return;
    }

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *const pHitParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            if (mcToTargetMCMap.find(pHitParticle) == mcToTargetMCMap.end())
                continue;

            MCParticleVector mcParticleContributionVector;
            for (const auto &mapEntry : pCaloHit->GetMCParticleWeightMap())
                mcParticleContributionVector.push_back(mapEntry.first);
            std::sort(mcParticleContributionVector.begin(), mcParticleContributionVector.end(), PointerLessThan<MCParticle>());

            MCParticleWeightMap targetWeightMap;
            for (const MCParticle *const pMCParticle : mcParticleContributionVector)
            {
                const float weight(pCaloHit->GetMCParticleWeightMap().at(pMCParticle));
                LArMCParticleHelper::MCRelationMap::const_iterator mcIter = mcToTargetMCMap.find(pMCParticle);

                if (mcToTargetMCMap.end() != mcIter)
                    targetWeightMap[mcIter->second] += weight;
            }

            MCParticleVector mcTargetContributionVector;
            for (const auto &mapEntry : targetWeightMap) mcTargetContributionVector.push_back(mapEntry.first);
            std::sort(mcTargetContributionVector.begin(), mcTargetContributionVector.end(), PointerLessThan<MCParticle>());

            float bestTargetWeight(0.f), targetWeightSum(0.f);

            for (const MCParticle *const pTargetMCParticle : mcTargetContributionVector)
            {
                const float targetWeight(targetWeightMap.at(pTargetMCParticle));
                targetWeightSum += targetWeight;

                if (targetWeight > bestTargetWeight)
                {
                    bestTargetWeight = targetWeight;
                }
            }

            if ((targetWeightSum < std::numeric_limits<float>::epsilon()) || ((bestTargetWeight / targetWeightSum) < minHitSharingFraction))
                continue;
            
            if(RejectBremsstrahlungHits(pCaloHit, leadingMCParticleToPostPhotonHitLists))
                continue;

            selectedCaloHitList.push_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArMuonLeadingHelper::RejectBremsstrahlungHits(const CaloHit *const pCaloHit, LeadingMCParticleToPostPhotonHitLists &leadingMCParticleToPostPhotonHitLists)
{
    const MCParticle *const pHitMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
    
    MCParticleList ancestorMCParticleList;
    LArMCParticleHelper::GetAllAncestorMCParticles(pHitMCParticle, ancestorMCParticleList);

    unsigned int highestTier(0);
    const MCParticle *leadingMCParticle(nullptr), *highestTierPhoton(nullptr);

    for (const MCParticle *const pAncestorMCParticle : ancestorMCParticleList)
    {
        if (LArMCParticleHelper::GetHierarchyTier(pAncestorMCParticle) == 1)
        {
            if (LArMuonLeadingHelper::IsLeading(pAncestorMCParticle))
                leadingMCParticle = pAncestorMCParticle;
        }

        if (pAncestorMCParticle->GetParticleId() == PHOTON)
        {
            if (LArMCParticleHelper::GetHierarchyTier(pAncestorMCParticle) > highestTier)
            {
                highestTier = LArMCParticleHelper::GetHierarchyTier(pAncestorMCParticle);
                highestTierPhoton = pAncestorMCParticle;
            }
        }
    }
   
    if (leadingMCParticle && highestTierPhoton)
        leadingMCParticleToPostPhotonHitLists[leadingMCParticle][highestTierPhoton].push_back(pCaloHit);

    if (leadingMCParticle && highestTierPhoton)
    {
        return true;
    }
    else
    {
        return false;
    }
}    

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::AddInReconstructablePostPhotonHits(const LeadingMCParticleToPostPhotonHitLists &leadingMCParticleToPostPhotonHitLists, const float maxBremsstrahlungSeparation,
                                                              LArMCParticleHelper::MCContributionMap &leadingMCToTrueHitListMap, const Pandora &/*pandora*/)
{
    MCParticleVector leadingMCParticleVector;
    for (auto &entry : leadingMCParticleToPostPhotonHitLists)
        leadingMCParticleVector.push_back(entry.first);
    std::sort(leadingMCParticleVector.begin(), leadingMCParticleVector.end(), LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pLeadingMCParticle : leadingMCParticleVector)
    {
        // DO NOT ADD IN HITS FOR WHICH THERE IS NO MAIN PARTICLE HITS
        if(leadingMCToTrueHitListMap.find(pLeadingMCParticle) == leadingMCToTrueHitListMap.end())
        {
            ////////////////////////////////
            /*
            std::cout << "ISOBEL - CASE WHERE THERE IS NO INITIAL DR TO ADD HITS ON TO" << std::endl;
            for (auto &entry : leadingMCParticleToPostPhotonHitLists.at(pLeadingMCParticle))
            {
                for (const CaloHit *const pCaloHit : entry.second)
                {
                    const CartesianVector &pos(pCaloHit->GetPositionVector());
                    PandoraMonitoringApi::AddMarkerToVisualization(pandora, &pos, "POSITION", BLACK, 2);
                }
            }
            */
            ////////////////////////////////            
            continue;
        }

        CaloHitList &leadingHitList(leadingMCToTrueHitListMap.at(pLeadingMCParticle));

        MCParticleVector pMCPhotonVector;
        for(auto &entry : leadingMCParticleToPostPhotonHitLists.at(pLeadingMCParticle))
            pMCPhotonVector.push_back(entry.first);
        std::sort(pMCPhotonVector.begin(), pMCPhotonVector.end(), LArMCParticleHelper::SortByMomentum);        

        for (const MCParticle *const pMCPhoton : pMCPhotonVector)
        {
            const CaloHitList postPhotonHitList(leadingMCParticleToPostPhotonHitLists.at(pLeadingMCParticle).at(pMCPhoton));

            if (LArClusterHelper::GetClosestDistanceWithShiftedHits(leadingHitList, postPhotonHitList) < maxBremsstrahlungSeparation)
                leadingHitList.insert(leadingHitList.begin(), postPhotonHitList.begin(), postPhotonHitList.end());
        }
    }

    //PandoraMonitoringApi::ViewEvent(pandora);
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
        const MCParticle *const pParentMCParticle(LArMCParticleHelper::GetParentMCParticle(pMCParticle));

        if (!LArMCParticleHelper::IsCosmicRay(pParentMCParticle))
            continue;

        if (pMCParticle == pParentMCParticle)
        {
            selectedParticles.push_back(pMCParticle);
        }
        else
        {
            if (LArMuonLeadingHelper::IsLeading(pMCParticle))
                selectedParticles.push_back(pMCParticle);
        }
    }

    std::sort(selectedParticles.begin(), selectedParticles.end(), LArMCParticleHelper::SortByMomentum);
}

    
} // namespace lar_content
