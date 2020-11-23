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
    {
        leadingMCParticleToPostPhotonHitLists[leadingMCParticle][pHitMCParticle].push_back(pCaloHit);        
        return true;
    }
    else
    {
        return false;
    }
}    

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::AddInReconstructablePostPhotonHits(const LeadingMCParticleToPostPhotonHitLists &leadingMCParticleToPostPhotonHitLists, const float maxBremsstrahlungSeparation,
    LArMCParticleHelper::MCContributionMap &leadingMCToTrueHitListMap, const Pandora &pandora)
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

        /*
        if (std::fabs(pLeadingMCParticle->GetEnergy() - 0.0296) < 0.01)
        {
        CaloHitList &leadingHitList(leadingMCToTrueHitListMap.at(pLeadingMCParticle));
        for(const CaloHit *const pCaloHit : leadingHitList)
        {
            const CartesianVector shiftedPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());
            PandoraMonitoringApi::AddMarkerToVisualization(pandora, &shiftedPosition, "CRL", RED, 2);
        }
        }
        */
            
        LArMuonLeadingHelper::AddHits(pLeadingMCParticle, leadingMCParticleToPostPhotonHitLists, maxBremsstrahlungSeparation, leadingMCToTrueHitListMap, TPC_VIEW_U, pandora);
        LArMuonLeadingHelper::AddHits(pLeadingMCParticle, leadingMCParticleToPostPhotonHitLists, maxBremsstrahlungSeparation, leadingMCToTrueHitListMap, TPC_VIEW_V, pandora);
        LArMuonLeadingHelper::AddHits(pLeadingMCParticle, leadingMCParticleToPostPhotonHitLists, maxBremsstrahlungSeparation, leadingMCToTrueHitListMap, TPC_VIEW_W, pandora);        
    }    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArMuonLeadingHelper::AddHits(const MCParticle *const pLeadingMCParticle, const LeadingMCParticleToPostPhotonHitLists &leadingMCParticleToPostPhotonHitLists,
                                   const float /*maxBremsstrahlungSeparation*/, LArMCParticleHelper::MCContributionMap &leadingMCToTrueHitListMap, const HitType &tpcView, const Pandora &pandora)
{
    CaloHitList leadingHitList;
    for (const CaloHit *const pCaloHit : leadingMCToTrueHitListMap.at(pLeadingMCParticle))
    {
        if (pCaloHit->GetHitType() == tpcView)
            leadingHitList.push_back(pCaloHit);
    }
    
    CaloHitList postBremsstrahlungHits;
    for (auto &entry : leadingMCParticleToPostPhotonHitLists.at(pLeadingMCParticle))
    {
        for (const CaloHit *const pCaloHit : entry.second)
        {
            if (pCaloHit->GetHitType() == tpcView)
                postBremsstrahlungHits.push_back(pCaloHit);
        }
    }

    //////////////////////////////
    for (const CaloHit *const pCaloHit : leadingHitList)
    {
        const CartesianVector shiftedPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());
        PandoraMonitoringApi::AddMarkerToVisualization(pandora, &shiftedPosition, "Main Delta Ray Hit", BLACK, 2);
    }

    PandoraMonitoringApi::ViewEvent(pandora);
    for (const CaloHit *const pCaloHit : postBremsstrahlungHits)
    {
        const CartesianVector shiftedPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());
        PandoraMonitoringApi::AddMarkerToVisualization(pandora, &shiftedPosition, "Post Bremsstrahlung Hit", BLUE, 2);
    }
        
    PandoraMonitoringApi::ViewEvent(pandora);
    //////////////////////////////    
    
    bool hitsAdded(true);
    while (hitsAdded)
    {
        hitsAdded = false;

        for (const CaloHit *const pPostBremsstrahlungHit : postBremsstrahlungHits)
        {
            if (std::find(leadingHitList.begin(), leadingHitList.end(), pPostBremsstrahlungHit) != leadingHitList.end())
                continue;
            
            const float separationDistance(LArClusterHelper::GetClosestDistanceWithShiftedHits(pPostBremsstrahlungHit, leadingHitList));

            if (separationDistance < 2.5f)
            {
                leadingHitList.push_back(pPostBremsstrahlungHit);
                hitsAdded = true;
                break;
            }
        }
    }

    for (const CaloHit *const pCaloHit : leadingHitList)
    {
        const CartesianVector shiftedPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());
        PandoraMonitoringApi::AddMarkerToVisualization(pandora, &shiftedPosition, "Reconstructable Delta Ray Hits", RED, 2);
    }

    PandoraMonitoringApi::ViewEvent(pandora);
}


/*
void LArMuonLeadingHelper::AddHits(const MCParticle *const pLeadingMCParticle, const LeadingMCParticleToPostPhotonHitLists &leadingMCParticleToPostPhotonHitLists,
    const float maxBremsstrahlungSeparation, LArMCParticleHelper::MCContributionMap &leadingMCToTrueHitListMap, const HitType &tpcView, const Pandora &pandora)
{
       MCParticleVector pMCPhotonVector;
        for(auto &entry : leadingMCParticleToPostPhotonHitLists.at(pLeadingMCParticle))
            pMCPhotonVector.push_back(entry.first);
        std::sort(pMCPhotonVector.begin(), pMCPhotonVector.end(), LArMCParticleHelper::SortByMomentum);
        
        bool hitsAdded(true);
        CaloHitList &leadingHitList(leadingMCToTrueHitListMap.at(pLeadingMCParticle));
        while (hitsAdded)
        {
            hitsAdded = false;
            
            CaloHitList leadingHits;
            for (const CaloHit *const star : leadingHitList)
            {
                if (star->GetHitType() == tpcView)
                    leadingHits.push_back(star);
            }

            CaloHitList hitsToAdd;
            for (const MCParticle *const pMCPhoton : pMCPhotonVector)
            {
                const CaloHitList postPhotonHitList(leadingMCParticleToPostPhotonHitLists.at(pLeadingMCParticle).at(pMCPhoton));

                CaloHitList postPhotonHits;
                for (const CaloHit *const star : postPhotonHitList)
                {
                    if (star->GetHitType() == tpcView)
                        postPhotonHits.push_back(star);
                }

                const float separationDistance(LArClusterHelper::GetClosestDistanceWithShiftedHits(leadingHits, postPhotonHits));

                /////////////////////////////////////
                for (const CaloHit *const pCaloHit : leadingHits)
                {
                    const CartesianVector shiftedPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());
                    PandoraMonitoringApi::AddMarkerToVisualization(pandora, &shiftedPosition, "Main Delta Ray Hit", RED, 2);
                }
                
                for (const CaloHit *const pCaloHit : postPhotonHits)
                {
                    const CartesianVector shiftedPosition(pCaloHit->GetPositionVector().GetX() - pCaloHit->GetX0(), pCaloHit->GetPositionVector().GetY(), pCaloHit->GetPositionVector().GetZ());
                    PandoraMonitoringApi::AddMarkerToVisualization(pandora, &shiftedPosition, "Post Photon Hit", BLUE, 2);
                }
                std::cout << "Separation Distance: " << separationDistance << std::endl;
                PandoraMonitoringApi::ViewEvent(pandora);
                /////////////////////////////////////
                
                if (separationDistance < maxBremsstrahlungSeparation)
                {
                    if ((separationDistance > 2.f) && (postPhotonHits.size() < 3))
                        continue;
                    
                    hitsToAdd.insert(hitsToAdd.begin(), postPhotonHits.begin(), postPhotonHits.end());
                    hitsAdded = true;
                    pMCPhotonVector.erase(std::find(pMCPhotonVector.begin(), pMCPhotonVector.end(), pMCPhoton));
                    break;
                }
            }

            leadingHitList.insert(leadingHitList.begin(), hitsToAdd.begin(), hitsToAdd.end());
        }
}
*/
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




    /*

        unsigned int highestTier(0);
        for (auto &entry : leadingMCParticleToPostPhotonHitLists.at(pLeadingMCParticle))
        {
            if (entry.first > highestTier)
                highestTier = entry.first;
        }

        for (unsigned int tier = 0; tier <= highestTier; ++tier)
        {
            if (leadingMCParticleToPostPhotonHitLists.at(pLeadingMCParticle).find(tier) == leadingMCParticleToPostPhotonHitLists.at(pLeadingMCParticle).end())
                continue;

            CaloHitList &leadingHitList(leadingMCToTrueHitListMap.at(pLeadingMCParticle));
            for (unsigned int isElectron = 0; isElectron < 2; ++isElectron)
            {
                CaloHitList leadingUHits, leadingVHits, leadingWHits;
                for (const CaloHit *const star : leadingHitList)
                {
                    if (star->GetHitType() == TPC_VIEW_U)
                        leadingUHits.push_back(star);

                    if (star->GetHitType() == TPC_VIEW_V)
                        leadingVHits.push_back(star);

                    if (star->GetHitType() == TPC_VIEW_W)
                        leadingWHits.push_back(star);
                }

                CaloHitList hitsToAdd;
                int pdg(isElectron ? 11 : -11);
                for (auto &entry : leadingMCParticleToPostPhotonHitLists.at(pLeadingMCParticle).at(tier))
                {
                    if (entry.first->GetParticleId() == pdg)
                    {
                        const CaloHitList postPhotonHitList(leadingMCParticleToPostPhotonHitLists.at(pLeadingMCParticle).at(tier).at(entry.first));
                        CaloHitList uHits, vHits, wHits;
                        for (const CaloHit *const star : postPhotonHitList)
                        {
                            if (star->GetHitType() == TPC_VIEW_U)
                                uHits.push_back(star);

                            if (star->GetHitType() == TPC_VIEW_V)
                                vHits.push_back(star);

                            if (star->GetHitType() == TPC_VIEW_W)
                                wHits.push_back(star);
                        }

                        if (LArClusterHelper::GetClosestDistanceWithShiftedHits(leadingUHits, uHits) < maxBremsstrahlungSeparation)
                            hitsToAdd.insert(hitsToAdd.begin(), uHits.begin(), uHits.end());

                        if (LArClusterHelper::GetClosestDistanceWithShiftedHits(leadingVHits, vHits) < maxBremsstrahlungSeparation)
                            hitsToAdd.insert(hitsToAdd.begin(), vHits.begin(), vHits.end());

                        if (LArClusterHelper::GetClosestDistanceWithShiftedHits(leadingWHits, wHits) < maxBremsstrahlungSeparation)
                            hitsToAdd.insert(hitsToAdd.begin(), wHits.begin(), wHits.end());
                    }
                }
                leadingHitList.insert(leadingHitList.begin(), hitsToAdd.begin(), hitsToAdd.end());
            }
        }
    }
}



          
                if (LArClusterHelper::GetClosestDistanceWithShiftedHits(leadingUHits, uHits) < maxBremsstrahlungSeparation)
                {
                if (std::fabs(pLeadingMCParticle->GetEnergy() - 0.3484) < 0.01)
                {
                PandoraMonitoringApi::VisualizeCaloHits(pandora, &leadingUHits, "LEADING", RED);
                PandoraMonitoringApi::VisualizeCaloHits(pandora, &uHits, "PHOTON", BLUE);                  
                std::cout << "U CLOSEST DISTANCE: " << LArClusterHelper::GetClosestDistanceWithShiftedHits(leadingUHits, uHits) << std::endl;
                PandoraMonitoringApi::ViewEvent(pandora);
                }
                hitsToAdd.insert(hitsToAdd.begin(), uHits.begin(), uHits.end());
            }

            if (LArClusterHelper::GetClosestDistanceWithShiftedHits(leadingVHits, vHits) < maxBremsstrahlungSeparation)
            {
                if (std::fabs(pLeadingMCParticle->GetEnergy() - 0.3484) < 0.01)
                {
                PandoraMonitoringApi::VisualizeCaloHits(pandora, &leadingVHits, "LEADING", RED);
                PandoraMonitoringApi::VisualizeCaloHits(pandora, &vHits, "PHOTON", BLUE);                  
                std::cout << "V CLOSEST DISTANCE: " << LArClusterHelper::GetClosestDistanceWithShiftedHits(leadingVHits, vHits) << std::endl;

                //for (const CaloHit *const sam : vHits)
                //std::cout << "HitParticle Vertex: " << MCParticleHelper::GetMainMCParticle(sam)->GetVertex() << std::endl;

                MCParticleList ancestorMCParticleList;
                LArMCParticleHelper::GetAllAncestorMCParticles(MCParticleHelper::GetMainMCParticle(vHits.front()), ancestorMCParticleList);
                for (const MCParticle *const pMCParticle : ancestorMCParticleList)
                {
                    std::cout << "ID: " << MCParticleHelper::GetMainMCParticle(vHits.front())->GetParticleId() << std::endl;
                    std::cout << "Tier: " << LArMCParticleHelper::GetHierarchyTier(pMCParticle) << ", " << pMCParticle->GetParticleId() << std::endl;
                }
                
                PandoraMonitoringApi::ViewEvent(pandora);
                }
                hitsToAdd.insert(hitsToAdd.begin(), vHits.begin(), vHits.end());
            }

            if (LArClusterHelper::GetClosestDistanceWithShiftedHits(leadingWHits, wHits) < maxBremsstrahlungSeparation)
            {
                if (std::fabs(pLeadingMCParticle->GetEnergy() - 0.3484) < 0.01)
                {
                PandoraMonitoringApi::VisualizeCaloHits(pandora, &leadingWHits, "LEADING", RED);
                PandoraMonitoringApi::VisualizeCaloHits(pandora, &wHits, "PHOTON", BLUE);                  
                std::cout << "W CLOSEST DISTANCE: " << LArClusterHelper::GetClosestDistanceWithShiftedHits(leadingWHits, wHits) << std::endl;
                PandoraMonitoringApi::ViewEvent(pandora);
                }
                hitsToAdd.insert(hitsToAdd.begin(), wHits.begin(), wHits.end());
            }                  
    */
