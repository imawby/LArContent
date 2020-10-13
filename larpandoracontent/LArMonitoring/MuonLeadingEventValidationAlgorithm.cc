/**
 *  @file   larpandoracontent/LArMonitoring/MuonLeadingEventValidationAlgorithm.cc
 *
 *  @brief  Implementation of the muon leading event validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArHelpers/LArMuonLeadingHelper.h"

#include "larpandoracontent/LArMonitoring/MuonLeadingEventValidationAlgorithm.h"

#include <sstream>

using namespace pandora;

namespace lar_content
{

MuonLeadingEventValidationAlgorithm::MuonLeadingEventValidationAlgorithm() :
    m_deltaRayMode(false),
    m_michelMode(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MuonLeadingEventValidationAlgorithm::~MuonLeadingEventValidationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::FillValidationInfo(const MCParticleList *const pMCParticleList, const CaloHitList *const pCaloHitList,
    const PfoList *const pPfoList, ValidationInfo &validationInfo) const
{
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);

    // First get unfolded matching (includes all particles within a CR hierarchy including the primary CR)
    this->SetUnfoldedMatching(pMCParticleList, pCaloHitList, pPfoList, validationInfo);

    if (pMCParticleList && pCaloHitList)
    {
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMuonLeadingHelper::SelectReconstructableLeadingParticles(pMCParticleList, pCaloHitList, m_validationParameters, targetMCParticleToHitsMap);

        LArMuonLeadingHelper::ValidationParameters validationParameters(m_validationParameters);
        validationParameters.m_minPrimaryGoodHits = 0;
        validationParameters.m_minHitsForGoodView = 0;
        validationParameters.m_minHitSharingFraction = 0.f;
        LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
        LArMuonLeadingHelper::SelectReconstructableLeadingParticles(pMCParticleList, pCaloHitList, validationParameters, allMCParticleToHitsMap);
        
        validationInfo.SetTargetMCParticleToHitsMap(targetMCParticleToHitsMap);
        validationInfo.SetAllMCParticleToHitsMap(allMCParticleToHitsMap);
    }
    
    if (pPfoList)
    {        
        PfoList nonMuonPfoList;
        LArMuonLeadingHelper::RemoveMuonPfosFromList(pPfoList, validationInfo.GetUnfoldedInterpretedMCToPfoHitSharingMap(), nonMuonPfoList);
        
        PfoList nonMuonLeadingPfoList;
        LArMuonLeadingHelper::SelectNonMuonLeadingPfos(nonMuonPfoList, nonMuonLeadingPfoList);
           
        LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
        LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(nonMuonLeadingPfoList, validationInfo.GetUnfoldedAllMCParticleToHitsMap(), pfoToHitsMap, true);

        validationInfo.SetPfoToHitsMap(pfoToHitsMap);
    }
    
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(validationInfo.GetPfoToHitsMap(), {validationInfo.GetAllMCParticleToHitsMap()}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);
    
    validationInfo.SetMCToPfoHitSharingMap(mcToPfoHitSharingMap);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMatching(validationInfo, interpretedMCToPfoHitSharingMap);
    validationInfo.SetInterpretedMCToPfoHitSharingMap(interpretedMCToPfoHitSharingMap);

    std::cout << "MC Hit Map size: " << validationInfo.GetTargetMCParticleToHitsMap().size() << std::endl;
    std::cout << "Pfo Hit Map size: " << validationInfo.GetPfoToHitsMap().size() << std::endl;
    std::cout << "Matching Map size: " << mcToPfoHitSharingMap.size() << std::endl;
    std::cout << "Interpretted Matching Map size: " << interpretedMCToPfoHitSharingMap.size() << std::endl;

}

//------------------------------------------------------------------------------------------------------------------------------------------
    
    void MuonLeadingEventValidationAlgorithm::ProcessOutput(const ValidationInfo &validationInfo, const bool /*useInterpretedMatching*/, const bool printToScreen, const bool fillTree) const
{
    // Unfolded hit ownership/sharing maps for cosmic ray muons and their descendents
    const LArMCParticleHelper::MCContributionMap &unfoldedTargetMCToHitsMap(validationInfo.GetUnfoldedTargetMCParticleToHitsMap());
    const LArMCParticleHelper::PfoContributionMap &unfoldedPfoToHitsMap(validationInfo.GetUnfoldedPfoToHitsMap());
    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &unfoldedMCToPfoHitSharingMap(validationInfo.GetUnfoldedInterpretedMCToPfoHitSharingMap());

    // Folded hit ownership/sharing maps for leading muon ionisation particles
    const LArMCParticleHelper::MCContributionMap &foldedAllMCToHitsMap(validationInfo.GetAllMCParticleToHitsMap());
    const LArMCParticleHelper::MCContributionMap &foldedTargetMCToHitsMap(validationInfo.GetTargetMCParticleToHitsMap());
    const LArMCParticleHelper::PfoContributionMap &foldedPfoToHitsMap(validationInfo.GetPfoToHitsMap());
    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &foldedMCToPfoHitSharingMap(validationInfo.GetInterpretedMCToPfoHitSharingMap());
    
    // Consider only delta rays from reconstructable CR muons
    MCParticleVector mcCRVector;
    for (auto &entry : unfoldedTargetMCToHitsMap)
    {
        if (LArMCParticleHelper::IsCosmicRay(entry.first))
            mcCRVector.push_back(entry.first);
    }
    std::sort(mcCRVector.begin(), mcCRVector.end(), LArMCParticleHelper::SortByMomentum);

    // Process matches
    int muonCount(0);    
    std::stringstream stringStream;
    for (const MCParticle *const pCosmicRay : mcCRVector)
    {
        // Cosmic ray parameters
        int muonId, nReconstructableChildCRLs(0);
        float mcE_CR, mcPX_CR, mcPY_CR, mcPZ_CR;
        int nMCHitsTotal_CR, nMCHitsU_CR, nMCHitsV_CR, nMCHitsW_CR;
        float mcVertexX_CR, mcVertexY_CR, mcVertexZ_CR, mcEndX_CR, mcEndY_CR, mcEndZ_CR;
        int nCorrectChildCRLs(0);
    
        // Leading particle parameters
        FloatVector mcE_CRL, mcPX_CRL, mcPY_CRL, mcPZ_CRL; 
        IntVector nMCHitsTotal_CRL, nMCHitsU_CRL, nMCHitsV_CRL, nMCHitsW_CRL;
        FloatVector mcVertexX_CRL, mcVertexY_CRL, mcVertexZ_CRL, mcEndX_CRL, mcEndY_CRL, mcEndZ_CRL;
        IntVector nAboveThresholdMatches_CRL, isCorrect_CRL, isCorrectParentLink_CRL;
        IntVector bestMatchNHitsTotal_CRL, bestMatchNHitsU_CRL, bestMatchNHitsV_CRL, bestMatchNHitsW_CRL;
        IntVector bestMatchNSharedHitsTotal_CRL, bestMatchNSharedHitsU_CRL, bestMatchNSharedHitsV_CRL, bestMatchNSharedHitsW_CRL;
        IntVector bestMatchNParentTrackHitsTotal_CRL, bestMatchNParentTrackHitsU_CRL, bestMatchNParentTrackHitsV_CRL, bestMatchNParentTrackHitsW_CRL;
        IntVector bestMatchNOtherTrackHitsTotal_CRL, bestMatchNOtherTrackHitsU_CRL, bestMatchNOtherTrackHitsV_CRL, bestMatchNOtherTrackHitsW_CRL;
        IntVector bestMatchNOtherShowerHitsTotal_CRL, bestMatchNOtherShowerHitsU_CRL, bestMatchNOtherShowerHitsV_CRL, bestMatchNOtherShowerHitsW_CRL;
        IntVector totalCRLHitsInBestMatchParentCR_CRL, uCRLHitsInBestMatchParentCR_CRL, vCRLHitsInBestMatchParentCR_CRL, wCRLHitsInBestMatchParentCR_CRL;

        // Move on if cosmic ray has not been reconstructed
        if (unfoldedMCToPfoHitSharingMap.at(pCosmicRay).empty())
            continue;

        // Obtain reconstructable leading particles
        MCParticleVector childLeadingParticles;
        for (const MCParticle *const pMuonChild : pCosmicRay->GetDaughterList())
        {
            if (m_deltaRayMode)
            {
                if (!LArMuonLeadingHelper::IsDeltaRay(pMuonChild))
                    continue;
            }

            if (m_michelMode)
            {
                if (!LArMuonLeadingHelper::IsMichel(pMuonChild))
                    continue;
            }

            if (!(m_michelMode || m_deltaRayMode))
                continue;

            // Move on if leading particle is not reconstructable
            LArMCParticleHelper::MCContributionMap::const_iterator iter(foldedTargetMCToHitsMap.find(pMuonChild));
            if (iter == foldedTargetMCToHitsMap.end())
                continue;

            childLeadingParticles.push_back(pMuonChild);
        }

        // Move on if cosmic ray has no leading delta ray child particles
        if (childLeadingParticles.empty())
            continue;
        
        std::sort(childLeadingParticles.begin(), childLeadingParticles.end(), LArMCParticleHelper::SortByMomentum);

        ++mounCount;

        // Pull cosmic ray info
        const CaloHitList &cosmicRayHitList(validationInfo.GetUnfoldedAllMCParticleToHitsMap().at(pCosmicRay));
        
        std::cout << "MC COSMIC RAY HITS" << std::endl;
        this->PrintHits(cosmicRayHitList, "MC CR", BLUE);

        muonId = muonCount;
        mcE_CR = pCosmicRay->GetEnergy();
        mcPX_CR = pCosmicRay->GetMomentum().GetX();
        mcPY_CR = pCosmicRay->GetMomentum().GetY();
        mcPZ_CR = pCosmicRay->GetMomentum().GetZ();        
        mcVertexX_CR = pCosmicRay->GetVertex().GetX();
        mcVertexY_CR = pCosmicRay->GetVertex().GetY();
        mcVertexZ_CR = pCosmicRay->GetVertex().GetZ();                
        mcEndX_CR = pCosmicRay->GetEndpoint().GetX();
        mcEndY_CR = pCosmicRay->GetEndpoint().GetY();
        mcEndZ_CR = pCosmicRay->GetEndpoint().GetZ();
        nMCHitsTotal_CR = cosmicRayHitList.size();
        nMCHitsU_CR = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, cosmicRayHitList);
        nMCHitsV_CR = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, cosmicRayHitList);        
        nMCHitsW_CR = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, cosmicRayHitList);
        nReconstructableChildCRLs = childLeadingParticles.size();
        
        stringStream << "\033[34m" << "(Parent CR: " << muonCount << ") " << "\033[0m" 
            << "Energy " << pCosmicRay->GetEnergy()
            << ", Dist. " << (pCosmicRay->GetEndpoint() - pCosmicRay->GetVertex()).GetMagnitude()
            << ", nMCHits " << cosmicRayHitList.size()
            << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, cosmicRayHitList)
            << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, cosmicRayHitList)
            << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, cosmicRayHitList) << ")"
            << ", nReconstructableCRLs " << nReconstructableChildCRLs << std::endl;   
        
        // Pull delta ray data 
        for (const MCParticle *const pLeadingParticle : childLeadingParticles)
        {
            // Pull delta ray MC info
            const CaloHitList &leadingParticleHitList(foldedAllMCToHitsMap.at(pLeadingParticle));

            std::cout << "MC DELTA RAY HITS" << std::endl;
            this->PrintHits(leadingParticleHitList, "MC DR", RED);           
            
            mcE_CRL.push_back(pLeadingParticle->GetEnergy());
            mcPX_CRL.push_back(pLeadingParticle->GetMomentum().GetX());
            mcPY_CRL.push_back(pLeadingParticle->GetMomentum().GetY());
            mcPZ_CRL.push_back(pLeadingParticle->GetMomentum().GetZ());
            mcVertexX_CRL.push_back(pLeadingParticle->GetVertex().GetX());
            mcVertexY_CRL.push_back(pLeadingParticle->GetVertex().GetY());
            mcVertexZ_CRL.push_back(pLeadingParticle->GetVertex().GetZ());        
            mcEndX_CRL.push_back(pLeadingParticle->GetEndpoint().GetX());
            mcEndY_CRL.push_back(pLeadingParticle->GetEndpoint().GetY());
            mcEndZ_CRL.push_back(pLeadingParticle->GetEndpoint().GetZ());   
            nMCHitsTotal_CRL.push_back(leadingParticleHitList.size());
            nMCHitsU_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, leadingParticleHitList));
            nMCHitsV_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, leadingParticleHitList));   
            nMCHitsW_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, leadingParticleHitList));

            stringStream << "\033[33m"  << "(Child "<< (m_deltaRayMode ? "DR)  " : "Michel)  ") << "\033[0m"
                << "Energy " << pLeadingParticle->GetEnergy()
                << ", Dist. " << (pLeadingParticle->GetEndpoint() - pLeadingParticle->GetVertex()).GetMagnitude()
                << ", nMCHits " << leadingParticleHitList.size()
                << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, leadingParticleHitList)
                << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, leadingParticleHitList)
                 << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, leadingParticleHitList) << ")" << std::endl;        

            // Look at the pfo matches
            int nMatches(0);
            int nAboveThresholdMatches(0);
            bool isBestMatchMatchedToCorrectCosmicRay(false);

            for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : foldedMCToPfoHitSharingMap.at(pLeadingParticle))
            {
                const ParticleFlowObject *const pMatchedPfo(pfoToSharedHits.first);
                const CaloHitList &pfoHitList(foldedPfoToHitsMap.at(pMatchedPfo));
                const CaloHitList &sharedHitList(pfoToSharedHits.second); 

                const bool isGoodMatch(this->IsGoodMatch(leadingParticleHitList, pfoHitList, sharedHitList));

                CaloHitList parentTrackHits, otherTrackHits, otherShowerHits;
                LArMuonLeadingHelper::GetPfoMatchContamination(pLeadingParticle, pfoHitList, parentTrackHits, otherTrackHits, otherShowerHits);

                // Check whether the reconstructed pfo has the correct parent-child link
                bool isMatchedToCorrectCosmicRay(false);
                const ParticleFlowObject *const pParentPfo(LArPfoHelper::GetParentPfo(pMatchedPfo));

                for (LArMCParticleHelper::PfoToSharedHitsVector::const_iterator cosmicRayMatchedPfoPair = unfoldedMCToPfoHitSharingMap.at(pCosmicRay).begin();
                    cosmicRayMatchedPfoPair != unfoldedMCToPfoHitSharingMap.at(pCosmicRay).end(); ++cosmicRayMatchedPfoPair)
                {
                    if (cosmicRayMatchedPfoPair->first == pParentPfo)
                    {
                        isMatchedToCorrectCosmicRay = true;
                        break;
                    }
                }                    

                CaloHitList leadingParticleHitsInParentCosmicRay;
                if (pParentPfo != pMatchedPfo)
                    LArMuonLeadingHelper::GetMuonPfoContaminationContribution(unfoldedPfoToHitsMap.at(pParentPfo), foldedAllMCToHitsMap.at(pLeadingParticle), leadingParticleHitsInParentCosmicRay);
                
                if (0 == nMatches++)
                {
                    isBestMatchMatchedToCorrectCosmicRay = isMatchedToCorrectCosmicRay;

                    if (isBestMatchMatchedToCorrectCosmicRay)
                    {
                        isCorrectParentLink_CRL.push_back(1);
                    }
                    else
                    {
                        isCorrectParentLink_CRL.push_back(0);
                    }                    
                    
                    bestMatchNHitsTotal_CRL.push_back(pfoHitList.size());
                    bestMatchNHitsU_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList));
                    bestMatchNHitsV_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList));
                    bestMatchNHitsW_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList));
                    
                    bestMatchNSharedHitsTotal_CRL.push_back(sharedHitList.size());
                    bestMatchNSharedHitsU_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList));
                    bestMatchNSharedHitsV_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList));
                    bestMatchNSharedHitsW_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList));

                    bestMatchNParentTrackHitsTotal_CRL.push_back(parentTrackHits.size());
                    bestMatchNParentTrackHitsU_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, parentTrackHits));
                    bestMatchNParentTrackHitsV_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, parentTrackHits));
                    bestMatchNParentTrackHitsW_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, parentTrackHits));

                    bestMatchNOtherTrackHitsTotal_CRL.push_back(otherTrackHits.size());
                    bestMatchNOtherTrackHitsU_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, otherTrackHits));
                    bestMatchNOtherTrackHitsV_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, otherTrackHits));
                    bestMatchNOtherTrackHitsW_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, otherTrackHits));

                    bestMatchNOtherShowerHitsTotal_CRL.push_back(otherShowerHits.size());
                    bestMatchNOtherShowerHitsU_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, otherShowerHits));
                    bestMatchNOtherShowerHitsV_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, otherShowerHits));
                    bestMatchNOtherShowerHitsW_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, otherShowerHits));
                    
                    totalCRLHitsInBestMatchParentCR_CRL.push_back(leadingParticleHitsInParentCosmicRay.size());
                    uCRLHitsInBestMatchParentCR_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, leadingParticleHitsInParentCosmicRay));
                    vCRLHitsInBestMatchParentCR_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, leadingParticleHitsInParentCosmicRay));
                    wCRLHitsInBestMatchParentCR_CRL.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, leadingParticleHitsInParentCosmicRay));
                }

                if (isGoodMatch)
                    ++nAboveThresholdMatches;
                                      
                stringStream << "-" << (!isGoodMatch ? "(Below threshold) " : "")
                    << "nPfoHits " << pfoHitList.size()
                    << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList) << ")"
                    << ", nMatchedHits " << sharedHitList.size()
                    << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList) << ")"
                    << ",  nCRLHitsInParentCR " << leadingParticleHitsInParentCosmicRay.size()
                    << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, leadingParticleHitsInParentCosmicRay)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, leadingParticleHitsInParentCosmicRay)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, leadingParticleHitsInParentCosmicRay) << ")" << std::endl                 
                    <<  (!isGoodMatch ? "                   " : " ")
                    << "nParentTrackHits " << parentTrackHits.size()
                    << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, parentTrackHits)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, parentTrackHits)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, parentTrackHits) << ")"
                    << ", nOtherTrackHits " << otherTrackHits.size()
                    << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, otherTrackHits)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, otherTrackHits)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, otherTrackHits) << ")"
                    << ", nOtherShowerHits " << otherShowerHits.size()
                    << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, otherShowerHits)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, otherShowerHits)
                             << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, otherShowerHits) << ")" << std::endl                    
                    <<  (!isGoodMatch ? "                   " : " ")                    
                    << (isMatchedToCorrectCosmicRay ? "Correct" :  "Incorrect") << "\033[0m" << " parent link" << std::endl;

                std::cout << stringStream.str() << std::endl;
                std::cout << "DELTA RAY pfo HITS" << std::endl;
                this->PrintHits(pfoHitList, "DR PFO", VIOLET);

                if (pParentPfo != pMatchedPfo)
                {
                    std::cout << "PARENT PFO" << std::endl;
                    const CaloHitList &parentCRHits(unfoldedPfoToHitsMap.at(pParentPfo));
                    this->PrintHits(parentCRHits, "DR PARENT PFO", DARKGREEN);
                }
                
             }
                
            nAboveThresholdMatches_CRL.push_back(nAboveThresholdMatches);

            bool isCorrect((nAboveThresholdMatches == 1) && isBestMatchMatchedToCorrectCosmicRay);
            if (isCorrect)
            {
                ++nCorrectChildCRLs;
                isCorrect_CRL.push_back(1);
            }
            else
            {
                isCorrect_CRL.push_back(0);
            }

            if (foldedMCToPfoHitSharingMap.at(pLeadingParticle).empty())
            {
                stringStream << "-" << "No matched pfo" << std::endl;
                std::cout << stringStream.str() << std::endl;
                
                isCorrectParentLink_CRL.push_back(0);
                bestMatchNHitsTotal_CRL.push_back(0); bestMatchNHitsU_CRL.push_back(0); bestMatchNHitsV_CRL.push_back(0); bestMatchNHitsW_CRL.push_back(0);
                bestMatchNSharedHitsTotal_CRL.push_back(0); bestMatchNSharedHitsU_CRL.push_back(0); bestMatchNSharedHitsV_CRL.push_back(0); bestMatchNSharedHitsW_CRL.push_back(0);
                bestMatchNParentTrackHitsTotal_CRL.push_back(0); bestMatchNParentTrackHitsU_CRL.push_back(0); bestMatchNParentTrackHitsV_CRL.push_back(0); bestMatchNParentTrackHitsW_CRL.push_back(0);
                bestMatchNOtherTrackHitsTotal_CRL.push_back(0); bestMatchNOtherTrackHitsU_CRL.push_back(0), bestMatchNOtherTrackHitsV_CRL.push_back(0); bestMatchNOtherTrackHitsW_CRL.push_back(0);
                bestMatchNOtherShowerHitsTotal_CRL.push_back(0); bestMatchNOtherShowerHitsU_CRL.push_back(0); bestMatchNOtherShowerHitsV_CRL.push_back(0), bestMatchNOtherShowerHitsW_CRL.push_back(0);
                totalCRLHitsInBestMatchParentCR_CRL.push_back(0);
                uCRLHitsInBestMatchParentCR_CRL.push_back(0); vCRLHitsInBestMatchParentCR_CRL.push_back(0); wCRLHitsInBestMatchParentCR_CRL.push_back(0);
            }
            
            stringStream << nAboveThresholdMatches << " above threshold matches" << std::endl
                         << "Reconstruction is "
                         << (isCorrect ?  "\033[32m" : "\033[31m")
                         << (isCorrect ?  "CORRECT" : "INCORRECT") << "\033[0m" << std::endl;

            std::cout << stringStream.str() << std::endl;
        }

        if (fillTree)
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventNumber", m_eventNumber - 1));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcE_CR", mcE_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPX_CR", mcPX_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPY_CR", mcPY_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPZ_CR", mcPZ_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsTotal_CR", nMCHitsTotal_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsU_CR", nMCHitsU_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsV_CR", nMCHitsV_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsW_CR", nMCHitsW_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexX_CR", mcVertexX_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexY_CR", mcVertexY_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexZ_CR", mcVertexZ_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndX_CR", mcEndX_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndY_CR", mcEndY_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndZ_CR", mcEndZ_CR));            
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nReconstructableChildCRLs", nReconstructableChildCRLs));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nCorrectChildCRLs", nCorrectChildCRLs));

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcE_CRL", &mcE_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPX_CRL", &mcPX_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPY_CRL", &mcPY_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPZ_CRL", &mcPZ_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsTotal_CRL", &nMCHitsTotal_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsU_CRL", &nMCHitsU_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsV_CRL", &nMCHitsV_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsW_CRL", &nMCHitsW_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexX_CRL", &mcVertexX_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexY_CRL", &mcVertexY_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexZ_CRL", &mcVertexZ_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndX_CRL", &mcEndX_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndY_CRL", &mcEndY_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndZ_CRL", &mcEndZ_CRL));             
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nAboveThresholdMatches_CRL", &nAboveThresholdMatches_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrect_CRL", &isCorrect_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrectParentLink_CRL", &isCorrectParentLink_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNHitsTotal_CRL", &bestMatchNHitsTotal_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNHitsU_CRL", &bestMatchNHitsU_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNHitsV_CRL", &bestMatchNHitsV_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNHitsW_CRL", &bestMatchNHitsW_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNSharedHitsTotal_CRL", &bestMatchNSharedHitsTotal_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNSharedHitsU_CRL", &bestMatchNSharedHitsU_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNSharedHitsV_CRL", &bestMatchNSharedHitsV_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNSharedHitsW_CRL", &bestMatchNSharedHitsW_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNParentTrackHitsTotal_CRL", &bestMatchNParentTrackHitsTotal_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNParentTrackHitsU_CRL", &bestMatchNParentTrackHitsU_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNParentTrackHitsV_CRL", &bestMatchNParentTrackHitsV_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNParentTrackHitsW_CRL", &bestMatchNParentTrackHitsW_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNOtherTrackHitsTotal_CRL", &bestMatchNOtherTrackHitsTotal_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNOtherTrackHitsU_CRL", &bestMatchNOtherTrackHitsU_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNOtherTrackHitsV_CRL", &bestMatchNOtherTrackHitsV_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNOtherTrackHitsW_CRL", &bestMatchNOtherTrackHitsW_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNOtherShowerHitsTotal_CRL", &bestMatchNOtherShowerHitsTotal_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNOtherShowerHitsU_CRL", &bestMatchNOtherShowerHitsU_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNOtherShowerHitsV_CRL", &bestMatchNOtherShowerHitsV_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNOtherShowerHitsW_CRL", &bestMatchNOtherShowerHitsW_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "totalCRLHitsInBestMatchParentCR_CRL", &totalCRLHitsInBestMatchParentCR_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "uCRLHitsInBestMatchParentCR_CRL", &uCRLHitsInBestMatchParentCR_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "vCRLHitsInBestMatchParentCR_CRL", &vCRLHitsInBestMatchParentCR_CRL));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "wCRLHitsInBestMatchParentCR_CRL", &wCRLHitsInBestMatchParentCR_CRL));

            PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName.c_str()));
        }

        stringStream << "------------------------------------------------------------------------------------------------" << std::endl;      
        stringStream << nCorrectChildCRLs << " / " << nReconstructableChildCRLs << " CRLs correctly reconstructed" << std::endl;
        stringStream << "------------------------------------------------------------------------------------------------" << std::endl;
        stringStream << "------------------------------------------------------------------------------------------------" << std::endl;
    }
    
    if (printToScreen)
    {
        std::cout << "jam" << std::endl;
        std::cout << stringStream.str() << std::endl;
    }
    

}

//------------------------------------------------------------------------------------------------------------------------------------------    

void MuonLeadingEventValidationAlgorithm::SetUnfoldedMatching(const MCParticleList *pMCParticleList, const CaloHitList *pCaloHitList, const PfoList *pPfoList,
    ValidationInfo &validationInfo) const
{
    ValidationInfo unfoldedValidationInfo;
    if (pMCParticleList && pCaloHitList)
    {
        LArMCParticleHelper::PrimaryParameters unfoldedParameters(m_primaryParameters);
        unfoldedParameters.m_foldBackHierarchy = false;
        
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, unfoldedParameters, LArMCParticleHelper::IsCosmicRay, targetMCParticleToHitsMap);

        LArMCParticleHelper::PrimaryParameters parameters(unfoldedParameters);
        parameters.m_minPrimaryGoodHits = 0;
        parameters.m_minHitsForGoodView = 0;
        parameters.m_minHitSharingFraction = 0.f;
        
        LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, allMCParticleToHitsMap);
        
        validationInfo.SetUnfoldedTargetMCParticleToHitsMap(targetMCParticleToHitsMap);
        validationInfo.SetUnfoldedAllMCParticleToHitsMap(allMCParticleToHitsMap);
        unfoldedValidationInfo.SetTargetMCParticleToHitsMap(targetMCParticleToHitsMap);
        unfoldedValidationInfo.SetAllMCParticleToHitsMap(allMCParticleToHitsMap);
    }

    if (pPfoList)
    {
        LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
        LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(*pPfoList, unfoldedValidationInfo.GetAllMCParticleToHitsMap(), pfoToHitsMap, false);

        validationInfo.SetUnfoldedPfoToHitsMap(pfoToHitsMap);
        unfoldedValidationInfo.SetPfoToHitsMap(pfoToHitsMap);
    }

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(unfoldedValidationInfo.GetPfoToHitsMap(), {unfoldedValidationInfo.GetAllMCParticleToHitsMap()}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);

    unfoldedValidationInfo.SetMCToPfoHitSharingMap(mcToPfoHitSharingMap);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMatching(unfoldedValidationInfo, interpretedMCToPfoHitSharingMap);

    validationInfo.SetUnfoldedInterpretedMCToPfoHitSharingMap(interpretedMCToPfoHitSharingMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MuonLeadingEventValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinPrimaryGoodHits", m_validationParameters.m_minPrimaryGoodHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitsForGoodView", m_validationParameters.m_minHitsForGoodView));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinPrimaryGoodViews", m_validationParameters.m_minPrimaryGoodViews));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SelectInputHits", m_validationParameters.m_selectInputHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinHitSharingFraction", m_validationParameters.m_minHitSharingFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPhotonPropagation", m_validationParameters.m_maxPhotonPropagation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FoldToPrimaries", m_validationParameters.m_foldBackHierarchy));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DeltaRayMode", m_deltaRayMode));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MichelMode", m_michelMode));
    
    return EventValidationBaseAlgorithm::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::PrintHits(const CaloHitList caloHitList, const std::string &stringTag, const Color &colour) const
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        if (pCaloHit->GetHitType() == TPC_VIEW_U)
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, stringTag, colour, 2);
    }
    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        if (pCaloHit->GetHitType() == TPC_VIEW_V)
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, stringTag, colour, 2);
    }
    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        if (pCaloHit->GetHitType() == TPC_VIEW_W)
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, stringTag, colour, 2);
    }
    PandoraMonitoringApi::ViewEvent(this->GetPandora());

        /*
        CartesianVector endpoint1(pCosmicRay->GetVertex()), endpoint2(pCosmicRay->GetEndpoint());
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpoint1, "muon endpoint ", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpoint2, "muon endpoint 2", RED, 2);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */
}

} // namespace lar_content    


        ////////////////////////////////////////
        /*
        std::cout << "All MC Cosmic Ray Hits" << std::endl;
        for (const CaloHit *const pCaloHit : cosmicRayHitList)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            if (pCaloHit->GetHitType() == TPC_VIEW_U)
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "MC_CR_U", BLUE, 2);
        }
        PandoraMonitoringApi::ViewEvent(this->GetPandora());

        for (const CaloHit *const pCaloHit : cosmicRayHitList)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            if (pCaloHit->GetHitType() == TPC_VIEW_V)
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "MC_CR_V", BLUE, 2);
        }
        PandoraMonitoringApi::ViewEvent(this->GetPandora());

        for (const CaloHit *const pCaloHit : cosmicRayHitList)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            if (pCaloHit->GetHitType() == TPC_VIEW_W)
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "MC_CR_W", BLUE, 2);
        }
        PandoraMonitoringApi::ViewEvent(this->GetPandora());

        std::cout << "Particle Id: " << pCosmicRay->GetParticleId() << std::endl;

        CartesianVector endpoint1(pCosmicRay->GetVertex()), endpoint2(pCosmicRay->GetEndpoint());
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpoint1, "muon endpoint ", RED, 2);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpoint2, "muon endpoint 2", RED, 2);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
        */
        ////////////////////////////////////////

            ////////////////////////////////////////
            /*
            std::cout << "All MC Delta Ray Hits" << std::endl;
            for (const CaloHit *const pCaloHit : leadingParticleHitList)
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                if (pCaloHit->GetHitType() == TPC_VIEW_U)
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "MC_DR_U", RED, 2);
            }
            PandoraMonitoringApi::ViewEvent(this->GetPandora());

            for (const CaloHit *const pCaloHit : leadingParticleHitList)
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                if (pCaloHit->GetHitType() == TPC_VIEW_V)
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "MC_DR_V", RED, 2);
            }
            PandoraMonitoringApi::ViewEvent(this->GetPandora());

            for (const CaloHit *const pCaloHit : leadingParticleHitList)
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                if (pCaloHit->GetHitType() == TPC_VIEW_W)
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "MC_DR_W", RED, 2);
            }
            PandoraMonitoringApi::ViewEvent(this->GetPandora());

            std::cout << "Particle Id: " << pCosmicRay->GetParticleId() << std::endl;

            CartesianVector endpoint1D(pLeadingParticle->GetVertex()), endpoint2D(pLeadingParticle->GetEndpoint());
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpoint1D, "delta endpoint ", VIOLET, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpoint2D, "delta endpoint 2", VIOLET, 2);
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
            */
            ////////////////////////////////////////


                ////////////////////////////////////////
                /*
                std::cout << "Pfo Delta Ray Hits" << std::endl;
                for (const CaloHit *const pCaloHit : pfoHitList)
                {
                    const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                    if (pCaloHit->GetHitType() == TPC_VIEW_U)
                        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "PFO_DR_U", VIOLET, 2);
                }
                PandoraMonitoringApi::ViewEvent(this->GetPandora());

                for (const CaloHit *const pCaloHit : pfoHitList)
                {
                    const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                    if (pCaloHit->GetHitType() == TPC_VIEW_V)
                        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "PFO_DR_V", VIOLET, 2);
                }
                PandoraMonitoringApi::ViewEvent(this->GetPandora());

                for (const CaloHit *const pCaloHit : pfoHitList)
                {
                    const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                    if (pCaloHit->GetHitType() == TPC_VIEW_W)
                        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "PFO_DR_W", VIOLET, 2);
                }
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
                */
                ////////////////////////////////////////

                ////////////////////////////////////////
                /*
                if (pParentPfo != pMatchedPfo)
                {
                    std::cout << "Pfo Parent Cosmic Ray Hits" << std::endl;
                    const CaloHitList &parentCRHits(unfoldedPfoToHitsMap.at(pParentPfo));
                    
                    for (const CaloHit *const pCaloHit : parentCRHits)
                    {   
                        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                        if (pCaloHit->GetHitType() == TPC_VIEW_U)
                            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "PFO_CR_U", DARKGREEN, 2);
                    }
                    PandoraMonitoringApi::ViewEvent(this->GetPandora());

                    for (const CaloHit *const pCaloHit : parentCRHits)
                    {
                        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                        if (pCaloHit->GetHitType() == TPC_VIEW_V)
                            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "PFO_CR_V", DARKGREEN, 2);
                    }
                    PandoraMonitoringApi::ViewEvent(this->GetPandora());

                    for (const CaloHit *const pCaloHit : parentCRHits)
                    {
                        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                        if (pCaloHit->GetHitType() == TPC_VIEW_W)
                            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "PFO_CR_W", DARKGREEN, 2);
                    }
                    PandoraMonitoringApi::ViewEvent(this->GetPandora());
                }
                else
                {
                    std::cout << "Pfo does not have a parent pfo" << std::endl;
                }
                */
                ////////////////////////////////////////   
