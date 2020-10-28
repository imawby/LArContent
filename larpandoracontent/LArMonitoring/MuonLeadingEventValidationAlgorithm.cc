/**
 *  @file   larpandoracontent/LArMonitoring/MuonLeadingEventValidationAlgorithm.cc
 *
 *  @brief  Implementation of the muon leading event validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
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
    m_michelMode(false),
    m_muonsToSkip(0),
    m_visualize(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

MuonLeadingEventValidationAlgorithm::~MuonLeadingEventValidationAlgorithm()
{
    if (m_writeToTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "OtherShowerContamination", m_fileName.c_str(), "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "EventValidationBaseAlgorithm: Unable to write tree " << "OtherShowerContamination" << " to file " << m_fileName << std::endl;
        }
    }
    if (m_writeToTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "OtherTrackContamination", m_fileName.c_str(), "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "EventValidationBaseAlgorithm: Unable to write tree " << "OtherTrackContamination" << " to file " << m_fileName << std::endl;
        }
    }
    if (m_writeToTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ParentTrackContamination", m_fileName.c_str(), "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "EventValidationBaseAlgorithm: Unable to write tree " << "ParentTrackContamination" << " to file " << m_fileName << std::endl;
        }
    }
    if (m_writeToTree)
    {
        try
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "CRLHitsInCR", m_fileName.c_str(), "UPDATE"));
        }
        catch (const StatusCodeException &)
        {
            std::cout << "EventValidationBaseAlgorithm: Unable to write tree " << "CRLHitsInCR" << " to file " << m_fileName << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::FillValidationInfo(const MCParticleList *const pMCParticleList, const CaloHitList *const pCaloHitList,
    const PfoList *const pPfoList, ValidationInfo &validationInfo) const
{
    //PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);
    //this->Investigate(pMCParticleList, pCaloHitList, *pPfoList);

    if (pMCParticleList && pCaloHitList)
    {
        // Get reconstructable MCParticle hit ownership map (non-muon leading hierarchy is folded whilst muon is unfolded)
        LArMuonLeadingHelper::ValidationParameters recoValidationParams(m_validationParameters);
        recoValidationParams.m_minHitSharingFraction = 0.7f;//0.90f;
        recoValidationParams.m_maxBremsstrahlungSeparation = 5.f; //0.f;//5.f;
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;        
        LArMuonLeadingHelper::SelectReconstructableLeadingParticles(pMCParticleList, pCaloHitList, recoValidationParams, targetMCParticleToHitsMap, this->GetPandora());    
        //validationInfo.SetTargetMCParticleToHitsMap(targetMCParticleToHitsMap);
        
        LArMuonLeadingHelper::ValidationParameters allValidationParams(m_validationParameters);
        allValidationParams.m_minPrimaryGoodHits = 0;
        allValidationParams.m_minHitsForGoodView = 0;
        allValidationParams.m_minHitSharingFraction = 0.7f;//0.90f;
        allValidationParams.m_maxBremsstrahlungSeparation = 5.f; //0.f;//5.f;
        LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
        LArMuonLeadingHelper::SelectReconstructableLeadingParticles(pMCParticleList, pCaloHitList, allValidationParams, allMCParticleToHitsMap, this->GetPandora());
        
        validationInfo.SetTargetMCParticleToHitsMap(targetMCParticleToHitsMap);
        validationInfo.SetAllMCParticleToHitsMap(allMCParticleToHitsMap);
    }

    if (pPfoList)
    {        
        LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
        LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(*pPfoList, validationInfo.GetAllMCParticleToHitsMap(), pfoToHitsMap, false);

        validationInfo.SetPfoToHitsMap(pfoToHitsMap);
    }
    
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(validationInfo.GetPfoToHitsMap(), {validationInfo.GetAllMCParticleToHitsMap()}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);
    
    validationInfo.SetMCToPfoHitSharingMap(mcToPfoHitSharingMap);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMatching(validationInfo, interpretedMCToPfoHitSharingMap);
    validationInfo.SetInterpretedMCToPfoHitSharingMap(interpretedMCToPfoHitSharingMap);

    //std::cout << "MC Hit Map size: " << validationInfo.GetTargetMCParticleToHitsMap().size() << std::endl;
    //std::cout << "Pfo Hit Map size: " << validationInfo.GetPfoToHitsMap().size() << std::endl;
    //std::cout << "Matching Map size: " << mcToPfoHitSharingMap.size() << std::endl;
    //std::cout << "Interpretted Matching Map size: " << interpretedMCToPfoHitSharingMap.size() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void MuonLeadingEventValidationAlgorithm::ProcessOutput(const ValidationInfo &validationInfo, const bool /*useInterpretedMatching*/, const bool printToScreen, const bool fillTree) const
{
    // Unfolded hit ownership/sharing maps for cosmic ray muons and their descendents
    //const LArMCParticleHelper::MCContributionMap &unfoldedTargetMCToHitsMap(validationInfo.GetUnfoldedTargetMCParticleToHitsMap());
    //const LArMCParticleHelper::PfoContributionMap &unfoldedPfoToHitsMap(validationInfo.GetUnfoldedPfoToHitsMap());
    //const LArMCParticleHelper::MCParticleToPfoHitSharingMap &unfoldedMCToPfoHitSharingMap(validationInfo.GetUnfoldedInterpretedMCToPfoHitSharingMap());

    // Folded hit ownership/sharing maps for leading muon ionisation particles
    const LArMCParticleHelper::MCContributionMap &foldedAllMCToHitsMap(validationInfo.GetAllMCParticleToHitsMap());
    const LArMCParticleHelper::MCContributionMap &foldedTargetMCToHitsMap(validationInfo.GetTargetMCParticleToHitsMap());
    const LArMCParticleHelper::PfoContributionMap &foldedPfoToHitsMap(validationInfo.GetPfoToHitsMap());
    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &foldedMCToPfoHitSharingMap(validationInfo.GetInterpretedMCToPfoHitSharingMap());
    
    // Consider only delta rays from reconstructable CR muons
    MCParticleVector mcCRVector;
    for (auto &entry : foldedTargetMCToHitsMap)
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
        int ID_CR, nReconstructableChildCRLs(0);
        float mcE_CR, mcPX_CR, mcPY_CR, mcPZ_CR;
        int nMCHitsTotal_CR, nMCHitsU_CR, nMCHitsV_CR, nMCHitsW_CR;
        float mcVertexX_CR, mcVertexY_CR, mcVertexZ_CR, mcEndX_CR, mcEndY_CR, mcEndZ_CR;
        int nCorrectChildCRLs(0);
    
        // Leading particle parameters
        FloatVector mcE_CRL, mcPX_CRL, mcPY_CRL, mcPZ_CRL;
        IntVector ID_CRL;        
        IntVector nMCHitsTotal_CRL, nMCHitsU_CRL, nMCHitsV_CRL, nMCHitsW_CRL;
        FloatVector mcVertexX_CRL, mcVertexY_CRL, mcVertexZ_CRL, mcEndX_CRL, mcEndY_CRL, mcEndZ_CRL;
        IntVector nAboveThresholdMatches_CRL, isCorrect_CRL, isCorrectParentLink_CRL, isBestMatchedCorrectParentLink_CRL;
        IntVector bestMatchNHitsTotal_CRL, bestMatchNHitsU_CRL, bestMatchNHitsV_CRL, bestMatchNHitsW_CRL;
        IntVector bestMatchNSharedHitsTotal_CRL, bestMatchNSharedHitsU_CRL, bestMatchNSharedHitsV_CRL, bestMatchNSharedHitsW_CRL;
        IntVector bestMatchNParentTrackHitsTotal_CRL, bestMatchNParentTrackHitsU_CRL, bestMatchNParentTrackHitsV_CRL, bestMatchNParentTrackHitsW_CRL;
        IntVector bestMatchNOtherTrackHitsTotal_CRL, bestMatchNOtherTrackHitsU_CRL, bestMatchNOtherTrackHitsV_CRL, bestMatchNOtherTrackHitsW_CRL;
        IntVector bestMatchNOtherShowerHitsTotal_CRL, bestMatchNOtherShowerHitsU_CRL, bestMatchNOtherShowerHitsV_CRL, bestMatchNOtherShowerHitsW_CRL;
        IntVector totalCRLHitsInBestMatchParentCR_CRL, uCRLHitsInBestMatchParentCR_CRL, vCRLHitsInBestMatchParentCR_CRL, wCRLHitsInBestMatchParentCR_CRL;

        // Contamination parameters
        IntVector  bestMatchOtherShowerHitsID_CRL, bestMatchOtherTrackHitsID_CRL, bestMatchParentTrackHitsID_CRL, bestMatchCRLHitsInCRID_CRL;
        FloatVector bestMatchOtherShowerHitsDistance_CRL, bestMatchOtherTrackHitsDistance_CRL, bestMatchParentTrackHitsDistance_CRL, bestMatchCRLHitsInCRDistance_CRL;

        // Move on if cosmic ray has not been reconstructed
        if (foldedMCToPfoHitSharingMap.at(pCosmicRay).empty())
            continue;

        // Obtain reconstructable leading particles
        MCParticleVector childLeadingParticles;
        for (const MCParticle *const pMuonChild : pCosmicRay->GetDaughterList())
        {
            if (!(m_michelMode || m_deltaRayMode))
                continue;
                        
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

        ++muonCount;

        if (muonCount < (m_muonsToSkip))
            continue;

        // Pull cosmic ray info
        const CaloHitList &cosmicRayHitList(foldedAllMCToHitsMap.at(pCosmicRay));

        ///////////////////////////////
        if(m_visualize)
        {
            std::cout << "MC COSMIC RAY HITS" << std::endl;
            this->PrintHits(cosmicRayHitList, "MC_CR", BLUE, false);
            CartesianVector endpoint1(pCosmicRay->GetVertex()), endpoint2(pCosmicRay->GetEndpoint());
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpoint1, "muon endpoint ", BLACK, 2);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpoint2, "muon endpoint 2", BLACK, 2);
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
        }
        ///////////////////////////////
        
        ID_CR = muonCount;
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
        int leadingCount(0);
        for (const MCParticle *const pLeadingParticle : childLeadingParticles)
        {
            // Pull delta ray MC info
            const CaloHitList &leadingParticleHitList(foldedAllMCToHitsMap.at(pLeadingParticle));
            ++leadingCount;

            ///////////////////////////////
            if (m_visualize)
            {
                std::cout << "MC DELTA RAY HITS" << std::endl;

                /*
                MCParticleList descendentMCParticleList;
                LArMCParticleHelper::GetAllDescendentMCParticles(pLeadingParticle, descendentMCParticleList);
                for (const MCParticle *const pDescendent : descendentMCParticleList)
                    std::cout << "MC Descendent: " << pDescendent->GetParticleId() << ", Tier: " << LArMCParticleHelper::GetHierarchyTier(pDescendent) << std::endl;
                std::cout << pLeadingParticle->GetEnergy() << std::endl;
                */
                this->PrintHits(leadingParticleHitList, "MC_DR", RED, false);
                CartesianVector endpoint1(pLeadingParticle->GetVertex()), endpoint2(pLeadingParticle->GetEndpoint());
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpoint1, "leading endpoint ", BLACK, 2);
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &endpoint2, "leading endpoint 2", BLACK, 2);
                PandoraMonitoringApi::ViewEvent(this->GetPandora());
            }
            ///////////////////////////////
            
            mcE_CRL.push_back(pLeadingParticle->GetEnergy());
            ID_CRL.push_back(leadingCount);
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

            stringStream << "\033[33m"  << "(Child "<< (m_deltaRayMode ? "DR: " : "Michel: ") << leadingCount << ")  " << "\033[0m"
                << "Energy " << pLeadingParticle->GetEnergy()
                << ", Dist. " << (pLeadingParticle->GetEndpoint() - pLeadingParticle->GetVertex()).GetMagnitude()
                << ", nMCHits " << leadingParticleHitList.size()
                << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, leadingParticleHitList)
                << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, leadingParticleHitList)
                 << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, leadingParticleHitList) << ")" << std::endl;        

            // Look at the pfo matches
            int nMatches(0), nAboveThresholdMatches(0);
            bool isCorrectParentLink(false);
            for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : foldedMCToPfoHitSharingMap.at(pLeadingParticle))
            {
                const ParticleFlowObject *const pMatchedPfo(pfoToSharedHits.first);
                const CaloHitList &pfoHitList(foldedPfoToHitsMap.at(pMatchedPfo));
                const CaloHitList &sharedHitList(pfoToSharedHits.second);

                const bool isGoodMatch(this->IsGoodMatch(leadingParticleHitList, pfoHitList, sharedHitList));

                ++nMatches;

                if (isGoodMatch)
                    ++nAboveThresholdMatches;

                CaloHitList parentTrackHits, otherTrackHits, otherShowerHits;
                LArMuonLeadingHelper::GetPfoMatchContamination(pLeadingParticle, pfoHitList, parentTrackHits, otherTrackHits, otherShowerHits);

                // Check whether the reconstructed pfo has the correct parent-child link
                bool isMatchedToCorrectCosmicRay(false);
                const ParticleFlowObject *const pParentPfo(LArPfoHelper::GetParentPfo(pMatchedPfo));

                for (LArMCParticleHelper::PfoToSharedHitsVector::const_iterator cosmicRayMatchedPfoPair = foldedMCToPfoHitSharingMap.at(pCosmicRay).begin();
                    cosmicRayMatchedPfoPair != foldedMCToPfoHitSharingMap.at(pCosmicRay).end(); ++cosmicRayMatchedPfoPair)
                {
                    if (cosmicRayMatchedPfoPair->first == pParentPfo)
                    {
                        isMatchedToCorrectCosmicRay = true;
                        break;
                    }
                }                    

                CaloHitList leadingParticleHitsInParentCosmicRay;
                if (isMatchedToCorrectCosmicRay)
                    LArMuonLeadingHelper::GetMuonPfoContaminationContribution(foldedPfoToHitsMap.at(pParentPfo), leadingParticleHitList, leadingParticleHitsInParentCosmicRay);
                
                if ((nAboveThresholdMatches == 1) && isGoodMatch)
                {
                    isCorrectParentLink = isMatchedToCorrectCosmicRay;
                        
                    if (isMatchedToCorrectCosmicRay)
                    {
                        isCorrectParentLink_CRL.push_back(1);
                    }
                    else
                    {
                        isCorrectParentLink_CRL.push_back(0);
                    }
                    
                    if (isMatchedToCorrectCosmicRay)
                    {
                        isBestMatchedCorrectParentLink_CRL.push_back(1);
                    }
                    else
                    {
                        isBestMatchedCorrectParentLink_CRL.push_back(0);
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

                    if (!otherShowerHits.empty())
                    {
                        bestMatchOtherShowerHitsID_CRL.insert(bestMatchOtherShowerHitsID_CRL.end(), otherShowerHits.size(), leadingCount);
                        this->FillContaminationHitsDistance(otherShowerHits, leadingParticleHitList, bestMatchOtherShowerHitsDistance_CRL);
                    }
                    if (!otherTrackHits.empty())
                    {
                        bestMatchOtherTrackHitsID_CRL.insert(bestMatchOtherTrackHitsID_CRL.end(), otherTrackHits.size(), leadingCount);
                        this->FillContaminationHitsDistance(otherTrackHits, leadingParticleHitList, bestMatchOtherTrackHitsDistance_CRL);
                    }
                    if (!parentTrackHits.empty())
                    {
                        bestMatchParentTrackHitsID_CRL.insert(bestMatchParentTrackHitsID_CRL.end(), parentTrackHits.size(), leadingCount);
                        this->FillContaminationHitsDistance(parentTrackHits, leadingParticleHitList, bestMatchParentTrackHitsDistance_CRL);
                    }

                    if (!leadingParticleHitsInParentCosmicRay.empty())
                    {
                        bestMatchCRLHitsInCRID_CRL.insert(bestMatchCRLHitsInCRID_CRL.end(), leadingParticleHitsInParentCosmicRay.size(), leadingCount);
                        this->FillContaminationHitsDistance(leadingParticleHitsInParentCosmicRay, cosmicRayHitList, bestMatchCRLHitsInCRDistance_CRL);
                    }
                }
                                      
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

                ///////////////////////////////
                if (m_visualize)
                {
                    std::cout << stringStream.str() << std::endl;
                    std::cout << "DELTA RAY PFO HITS" << std::endl;
                    this->PrintHits(pfoHitList, otherShowerHits, otherTrackHits, parentTrackHits, "DR_PFO");

                    if (pParentPfo != pMatchedPfo)
                    {
                        std::cout << "PARENT PFO" << std::endl;
                        const CaloHitList &parentCRHits(foldedPfoToHitsMap.at(pParentPfo));
                        this->PrintHits(parentCRHits, leadingParticleHitList, "DR_PARENT_PFO");
                    }
                }
                ///////////////////////////////
             }
                
            nAboveThresholdMatches_CRL.push_back(nAboveThresholdMatches);

            const bool isCorrect((nAboveThresholdMatches == 1) && isCorrectParentLink);
            
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

                if (m_visualize)
                {
                    std::cout << stringStream.str() << std::endl;
                }
            }

            if (nAboveThresholdMatches == 0)
            {
                isCorrectParentLink_CRL.push_back(0);                
                isBestMatchedCorrectParentLink_CRL.push_back(0);
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

            if (m_visualize)
                std::cout << stringStream.str() << std::endl;
        }

        if (fillTree)
        {
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventNumber", m_eventNumber - 1));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ID_CR", ID_CR));
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

            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "ID_CRL", &ID_CRL));
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
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isBestMatchedCorrectParentLink_CRL", &isBestMatchedCorrectParentLink_CRL));
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
            
            if (!bestMatchOtherShowerHitsID_CRL.empty())
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "OtherShowerContamination", "bestMatchOtherShowerHitsEventNumber", m_eventNumber - 1));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "OtherShowerContamination", "bestMatchOtherShowerHitsID_CR", ID_CR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "OtherShowerContamination", "bestMatchOtherShowerHitsID_CRL", &bestMatchOtherShowerHitsID_CRL));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "OtherShowerContamination", "bestMatchOtherShowerHitsDistance_CRL", &bestMatchOtherShowerHitsDistance_CRL));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), "OtherShowerContamination"));
            }

            if (!bestMatchOtherTrackHitsID_CRL.empty())
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "OtherTrackContamination", "bestMatchOtherTrackHitsEventNumber", m_eventNumber - 1));                
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "OtherTrackContamination", "bestMatchOtherTrackHitsID_CR", ID_CR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "OtherTrackContamination", "bestMatchOtherTrackHitsID_CRL", &bestMatchOtherTrackHitsID_CRL));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "OtherTrackContamination", "bestMatchOtherTrackHitsDistance_CRL", &bestMatchOtherTrackHitsDistance_CRL));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), "OtherTrackContamination"));
            }

            if (!bestMatchParentTrackHitsID_CRL.empty())
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ParentTrackContamination", "bestMatchParentTrackHitsEventNumber", m_eventNumber - 1));                
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ParentTrackContamination", "bestMatchParentTrackHitsID_CR", ID_CR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ParentTrackContamination", "bestMatchParentTrackHitsID_CRL", &bestMatchParentTrackHitsID_CRL));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "ParentTrackContamination", "bestMatchParentTrackHitsDistance_CRL", &bestMatchParentTrackHitsDistance_CRL));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), "ParentTrackContamination"));
            }

            if (!bestMatchCRLHitsInCRID_CRL.empty())
            {
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CRLHitsInCR", "bestMatchCRLHitsInCREventNumber", m_eventNumber -1));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CRLHitsInCR", "bestMatchCRLHitsInCRID_CR", ID_CR));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CRLHitsInCR", "bestMatchCRLHitsInCRID_CRL", &bestMatchCRLHitsInCRID_CRL));
                PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CRLHitsInCR", "bestMatchCRLHitsInCRDistance_CRL", &bestMatchCRLHitsInCRDistance_CRL));
                PANDORA_MONITORING_API(FillTree(this->GetPandora(), "CRLHitsInCR"));
            }
        }

        stringStream << "------------------------------------------------------------------------------------------------" << std::endl;      
        stringStream << nCorrectChildCRLs << " / " << nReconstructableChildCRLs << " CRLs correctly reconstructed" << std::endl;
        stringStream << "------------------------------------------------------------------------------------------------" << std::endl;
        stringStream << "------------------------------------------------------------------------------------------------" << std::endl;
    }
    
    if (printToScreen && !m_visualize)
    {
        std::cout << stringStream.str() << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------    
/*
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
*/
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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MuonsToSkip", m_muonsToSkip));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "Visualize", m_visualize));        
    
    return EventValidationBaseAlgorithm::ReadSettings(xmlHandle);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::PrintHits(const CaloHitList caloHitList, const std::string &stringTag, const Color &colour, bool print) const
{
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

        if (pCaloHit->GetHitType() == TPC_VIEW_U)
        {
            Color printColour(colour);
            if(print)
            {
                std::cout << "/////////////////////" << std::endl;                
                const MCParticleWeightMap &weightMap(pCaloHit->GetMCParticleWeightMap());
                for (auto &entry : weightMap)
                {
                    std::cout << "MCParticle: " << entry.first->GetParticleId() << ", Weight: " << entry.second << std::endl;
                    if ((std::abs(entry.first->GetParticleId()) == 13) && (entry.second > 0.3))
                        printColour = VIOLET;
                        
                }
            }
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, stringTag, printColour, 2);
            //PandoraMonitoringApi::Pause(this->GetPandora());
        }
    }
    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        
        if (pCaloHit->GetHitType() == TPC_VIEW_V)
        {
            Color printColour(colour);
            if(print)
            {
                std::cout << "/////////////////////" << std::endl;
                const MCParticleWeightMap &weightMap(pCaloHit->GetMCParticleWeightMap());
                for (auto &entry : weightMap)
                {
                    std::cout << "MCParticle: " << entry.first->GetParticleId() << ", Weight: " << entry.second << std::endl;
                    if ((std::abs(entry.first->GetParticleId()) == 13) && (entry.second > 0.3))
                        printColour = VIOLET;
                }
            }
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, stringTag, printColour, 2);
            //PandoraMonitoringApi::Pause(this->GetPandora());            
        }
    }
    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

        if (pCaloHit->GetHitType() == TPC_VIEW_W)
        {
            Color printColour(colour);            
            if(print)
            {
                std::cout << "/////////////////////" << std::endl;                
                const MCParticleWeightMap &weightMap(pCaloHit->GetMCParticleWeightMap());
                for (auto &entry : weightMap)
                {
                    std::cout << "MCParticle: " << entry.first->GetParticleId() << ", Weight: " << entry.second << std::endl;
                    if ((std::abs(entry.first->GetParticleId()) == 13) && (entry.second > 0.3))
                        printColour = VIOLET;                    
                }
            }
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, stringTag, printColour, 2);            
            //PandoraMonitoringApi::Pause(this->GetPandora());
        }
    }
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::PrintHits(const CaloHitList totalCaloHitList, const CaloHitList otherShowerCaloHitList,
    const CaloHitList otherTrackCaloHitList, const CaloHitList parentTrackCaloHitList, const std::string &stringTag) const
{
    for (const CaloHit *const pCaloHit : totalCaloHitList)
    {
        if (pCaloHit->GetHitType() != TPC_VIEW_U)
            continue;

        Color color(BLACK);
        std::string newStringTag(stringTag);
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        
        if (std::find(otherShowerCaloHitList.begin(), otherShowerCaloHitList.end(), pCaloHit) != otherShowerCaloHitList.end())
        {
            newStringTag += "_OTHER_SHOWER"; color = VIOLET;
        }
        if (std::find(otherTrackCaloHitList.begin(), otherTrackCaloHitList.end(), pCaloHit) != otherTrackCaloHitList.end())
        {
            newStringTag += "_OTHER_TRACK"; color = RED;
        }
        if (std::find(parentTrackCaloHitList.begin(), parentTrackCaloHitList.end(), pCaloHit) != parentTrackCaloHitList.end())
        {
            newStringTag += "_PARENT_TRACK"; color = BLUE;
        }
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, newStringTag, color, 2);
    }
    PandoraMonitoringApi::ViewEvent(this->GetPandora());  

    for (const CaloHit *const pCaloHit : totalCaloHitList)
    {
        if (pCaloHit->GetHitType() != TPC_VIEW_V)
            continue;

        Color color(BLACK);
        std::string newStringTag(stringTag);        
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            
        if (std::find(otherShowerCaloHitList.begin(), otherShowerCaloHitList.end(), pCaloHit) != otherShowerCaloHitList.end())
        {
            newStringTag += "_OTHER_SHOWER"; color = VIOLET;
        }
        if (std::find(otherTrackCaloHitList.begin(), otherTrackCaloHitList.end(), pCaloHit) != otherTrackCaloHitList.end())
        {
            newStringTag += "_OTHER_TRACK"; color = RED;
        }
        if (std::find(parentTrackCaloHitList.begin(), parentTrackCaloHitList.end(), pCaloHit) != parentTrackCaloHitList.end())
        {
            newStringTag += "_PARENT_TRACK"; color = BLUE;
        }
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, newStringTag, color, 2);
    }
    PandoraMonitoringApi::ViewEvent(this->GetPandora());  

    for (const CaloHit *const pCaloHit : totalCaloHitList)
    {
        if (pCaloHit->GetHitType() != TPC_VIEW_W)
            continue;

        Color color(BLACK);
        std::string newStringTag(stringTag);        
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            
        if (std::find(otherShowerCaloHitList.begin(), otherShowerCaloHitList.end(), pCaloHit) != otherShowerCaloHitList.end())
        {
            newStringTag += "_OTHER_SHOWER"; color = VIOLET;
        }
        if (std::find(otherTrackCaloHitList.begin(), otherTrackCaloHitList.end(), pCaloHit) != otherTrackCaloHitList.end())
        {
            newStringTag += "_OTHER_TRACK"; color = RED;
        }
        if (std::find(parentTrackCaloHitList.begin(), parentTrackCaloHitList.end(), pCaloHit) != parentTrackCaloHitList.end())
        {
            newStringTag += "_PARENT_TRACK"; color = BLUE;
        }
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, newStringTag, color, 2);
    }
    PandoraMonitoringApi::ViewEvent(this->GetPandora());  
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::PrintHits(const CaloHitList totalCaloHitList, const CaloHitList leadingCaloHitList,
    const std::string &stringTag) const
{
    for (const CaloHit *const pCaloHit : totalCaloHitList)
    {
        if (pCaloHit->GetHitType() != TPC_VIEW_U)
            continue;

        Color color(DARKGREEN);
        std::string newStringTag(stringTag);
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
        
        if (std::find(leadingCaloHitList.begin(), leadingCaloHitList.end(), pCaloHit) != leadingCaloHitList.end())
        {
            newStringTag += "_LEADING"; color = RED;
        }
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, newStringTag, color, 2);
    }
    PandoraMonitoringApi::ViewEvent(this->GetPandora());  

    for (const CaloHit *const pCaloHit : totalCaloHitList)
    {
        if (pCaloHit->GetHitType() != TPC_VIEW_V)
            continue;

        Color color(DARKGREEN);
        std::string newStringTag(stringTag);        
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            
        if (std::find(leadingCaloHitList.begin(), leadingCaloHitList.end(), pCaloHit) != leadingCaloHitList.end())
        {
            newStringTag += "_LEADING"; color = RED;
        }
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, newStringTag, color, 2);
    }
    PandoraMonitoringApi::ViewEvent(this->GetPandora());  

    for (const CaloHit *const pCaloHit : totalCaloHitList)
    {
        if (pCaloHit->GetHitType() != TPC_VIEW_W)
            continue;

        Color color(DARKGREEN);
        std::string newStringTag(stringTag);        
        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            
        if (std::find(leadingCaloHitList.begin(), leadingCaloHitList.end(), pCaloHit) != leadingCaloHitList.end())
        {
            newStringTag += "_LEADING"; color = RED;
        }
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, newStringTag, color, 2);
    }
    PandoraMonitoringApi::ViewEvent(this->GetPandora());  
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::FillContaminationHitsDistance(const CaloHitList &contaminationHits, const CaloHitList &leadingMCHits,
    FloatVector &bestMatchContaminationHitsDistance) const
{
    CaloHitList leadingU, leadingV, leadingW;
    this->GetHitsOfType(leadingMCHits, TPC_VIEW_U, leadingU);
    this->GetHitsOfType(leadingMCHits, TPC_VIEW_V, leadingV);
    this->GetHitsOfType(leadingMCHits, TPC_VIEW_W, leadingW);    
    
    for (const CaloHit *const pContaminationHit : contaminationHits)
    {
        const CartesianVector &hitPosition(pContaminationHit->GetPositionVector());
        
        if (pContaminationHit->GetHitType() == TPC_VIEW_U)
            bestMatchContaminationHitsDistance.push_back(LArClusterHelper::GetClosestDistance(hitPosition, leadingU));

        if (pContaminationHit->GetHitType() == TPC_VIEW_V)
            bestMatchContaminationHitsDistance.push_back(LArClusterHelper::GetClosestDistance(hitPosition, leadingV));
            
        if (pContaminationHit->GetHitType() == TPC_VIEW_W)
            bestMatchContaminationHitsDistance.push_back(LArClusterHelper::GetClosestDistance(hitPosition, leadingW));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::GetHitsOfType(const CaloHitList &inputList, const HitType &hitType, CaloHitList &outputList) const
{
    for (const CaloHit *const pCaloHit : inputList)
    {
        if (pCaloHit->GetHitType() == hitType)
            outputList.push_back(pCaloHit);
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void MuonLeadingEventValidationAlgorithm::Investigate(const MCParticleList *const pMCParticleList, const CaloHitList *const pCaloHitList, const PfoList &pfoList) const
{
    ////////////////////////////
    ValidationInfo validationInfo;
    
    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    LArMuonLeadingHelper::ValidationParameters jam(m_validationParameters);
    jam.m_minHitSharingFraction = 0.70f;
    jam.m_maxBremsstrahlungSeparation = 0.f;
    LArMuonLeadingHelper::SelectReconstructableLeadingParticles(pMCParticleList, pCaloHitList, jam, targetMCParticleToHitsMap, this->GetPandora());        
    validationInfo.SetTargetMCParticleToHitsMap(targetMCParticleToHitsMap);

    LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
    LArMuonLeadingHelper::ValidationParameters jamAll(m_validationParameters);
    jamAll.m_minPrimaryGoodHits = 0;
    jamAll.m_minHitsForGoodView = 0;    
    jamAll.m_minHitSharingFraction = 0.70f;
    jamAll.m_maxBremsstrahlungSeparation = std::numeric_limits<float>::max();
    LArMuonLeadingHelper::SelectReconstructableLeadingParticles(pMCParticleList, pCaloHitList, jamAll, allMCParticleToHitsMap, this->GetPandora());     
    validationInfo.SetAllMCParticleToHitsMap(allMCParticleToHitsMap);

    // To identify hits rejected because of their hit sharing fraction
    LArMuonLeadingHelper::ValidationParameters noHitFractionParameters(m_validationParameters);
    noHitFractionParameters.m_minPrimaryGoodHits = 0;
    noHitFractionParameters.m_minHitsForGoodView = 0;
    noHitFractionParameters.m_minHitSharingFraction = 0.f;
    noHitFractionParameters.m_maxBremsstrahlungSeparation = 0.f;        
    LArMCParticleHelper::MCContributionMap noHitFraction;
    LArMuonLeadingHelper::SelectReconstructableLeadingParticles(pMCParticleList, pCaloHitList, noHitFractionParameters, noHitFraction, this->GetPandora());

    // To identify hits rejected because of their bremsstrahlung propagation length
    LArMuonLeadingHelper::ValidationParameters noPhotonPropagationCut(m_validationParameters);
    noPhotonPropagationCut.m_minPrimaryGoodHits = 0;
    noPhotonPropagationCut.m_minHitsForGoodView = 0;
    noPhotonPropagationCut.m_minHitSharingFraction = 0.70f;
    noPhotonPropagationCut.m_maxBremsstrahlungSeparation = std::numeric_limits<float>::max();
    LArMCParticleHelper::MCContributionMap noPhotonPropagation;
    LArMuonLeadingHelper::SelectReconstructableLeadingParticles(pMCParticleList, pCaloHitList, noPhotonPropagationCut, noPhotonPropagation, this->GetPandora());    

    LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
    LArMCParticleHelper::MCContributionMapVector trueHitMaps = {noHitFraction, noPhotonPropagation}; //need to fill with hits would expect to get
    LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(pfoList, trueHitMaps, pfoToHitsMap, false);
    validationInfo.SetPfoToHitsMap(pfoToHitsMap);

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(pfoToHitsMap, {targetMCParticleToHitsMap}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);
    validationInfo.SetMCToPfoHitSharingMap(mcToPfoHitSharingMap);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMatching(validationInfo, interpretedMCToPfoHitSharingMap);
    validationInfo.SetInterpretedMCToPfoHitSharingMap(interpretedMCToPfoHitSharingMap);
    ////////////////////////////    




    //Print matches
        for (auto &entry : targetMCParticleToHitsMap)
        {
            if (std::fabs(entry.first->GetParticleId()) == 13)
            {
                if (targetMCParticleToHitsMap.find(entry.first) != targetMCParticleToHitsMap.end())
                {
                    const MCParticleList &childParticles(entry.first->GetDaughterList());
                    
                    for (const MCParticle *const pChild : childParticles)
                    {
                        if (targetMCParticleToHitsMap.find(pChild) != targetMCParticleToHitsMap.end())
                        {
                            
                            MCParticleList hierarchy;
                            LArMCParticleHelper::GetAllDescendentMCParticles(pChild, hierarchy);

                            for (const MCParticle *const pJam : hierarchy)
                                std::cout << "MCParticle: " << pJam->GetParticleId() << ", Hierarchy: " << LArMCParticleHelper::GetHierarchyTier(pJam) << std::endl;
                            
                            CaloHitList reconstructableHits;
                            for (const CaloHit *const pCaloHit : targetMCParticleToHitsMap.at(pChild))
                                reconstructableHits.push_back(pCaloHit);

                            CaloHitList sharingFractionRejectedHits;
                            if (noHitFraction.find(pChild) != noHitFraction.end())
                            {
                                for (const CaloHit *const pCaloHit : noHitFraction.at(pChild))
                                {
                                    if (std::find(reconstructableHits.begin(), reconstructableHits.end(), pCaloHit) == reconstructableHits.end())
                                        sharingFractionRejectedHits.push_back(pCaloHit);
                                }
                            }

                            CaloHitList photonPropagationRejectedHits;
                            if (noPhotonPropagation.find(pChild) != noPhotonPropagation.end())
                            {
                                for (const CaloHit *const pCaloHit : noPhotonPropagation.at(pChild))
                                {
                                    if (std::find(reconstructableHits.begin(), reconstructableHits.end(), pCaloHit) == reconstructableHits.end())
                                        photonPropagationRejectedHits.push_back(pCaloHit);
                                }
                            }                            

                            for (const CaloHit *const pCaloHit : reconstructableHits)
                            {
                                if (pCaloHit->GetHitType() != TPC_VIEW_U)
                                    continue;
                                
                                const CartesianVector &position(pCaloHit->GetPositionVector());
                                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "RECO", RED, 2);
                            }

                            for (const CaloHit *const pCaloHit : sharingFractionRejectedHits)
                            {
                                if (pCaloHit->GetHitType() != TPC_VIEW_U)
                                    continue;

                                std::cout << "/////////////////////" << std::endl;
                                const MCParticleWeightMap &weightMap(pCaloHit->GetMCParticleWeightMap());
                                for (auto &frog : weightMap)
                                {
                                    std::cout << "MCParticle: " << frog.first->GetParticleId() << ", Weight: " << frog.second << std::endl;
                                }

                                const CartesianVector &position(pCaloHit->GetPositionVector());
                                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "SHARING FRACTION", VIOLET, 2);
                            }

                            for (const CaloHit *const pCaloHit : photonPropagationRejectedHits)
                            {
                                if (pCaloHit->GetHitType() != TPC_VIEW_U)
                                    continue;
                                
                                const CartesianVector &position(pCaloHit->GetPositionVector());
                                //const CartesianVector translatedPosition(position.GetX() - pCaloHit->GetX0(), position.GetY(), position.GetZ());

                                std::cout << "Distance: " << LArClusterHelper::GetClosestDistanceWithShiftedHits(pCaloHit, targetMCParticleToHitsMap.at(pChild)) << std::endl;
                                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "PHOTON LENGTH", BLUE, 2);
                                //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &translatedPosition, "PHOTON TRANSLATED", DARKBLUE, 2);                                
                                //std::cout << "X0: " << pCaloHit->GetX0() << std::endl;
                            }                                               

                            PandoraMonitoringApi::ViewEvent(this->GetPandora());
                        
                            for (auto &pfoMatch : interpretedMCToPfoHitSharingMap.at(pChild))
                            {
                                for (const CaloHit *const pfoHit : pfoToHitsMap.at(pfoMatch.first))
                                {
                                    const CartesianVector &position(pfoHit->GetPositionVector());
                                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "PFO", DARKGREEN, 2);
                                }

                                PandoraMonitoringApi::ViewEvent(this->GetPandora());
                            }
                        }                            
                    }
                }
            }
        }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool MuonLeadingEventValidationAlgorithm::PrintPhoton(const MCParticle *const pOriginalPrimary, const MCParticle *const pThisMCParticle,
    const MCParticle *const pHitMCParticle, const float maxPhotonPropagation, const HitType &hitType) const
{
    if (NEUTRON == std::abs(pThisMCParticle->GetParticleId()))
        return false;

    MCParticleList hierarchy;
    LArMCParticleHelper::GetAllAncestorMCParticles(pHitMCParticle, hierarchy);

    if ((PHOTON == pThisMCParticle->GetParticleId()) && (PHOTON != pOriginalPrimary->GetParticleId()) && (E_MINUS != std::abs(pOriginalPrimary->GetParticleId())))
    {
        const CartesianVector &vertex(pThisMCParticle->GetVertex()), &endpoint(pThisMCParticle->GetEndpoint());
        const CartesianVector projectedVertex(LArGeometryHelper::ProjectPosition(this->GetPandora(), vertex, hitType));
        const CartesianVector projectedEndpoint(LArGeometryHelper::ProjectPosition(this->GetPandora(), endpoint, hitType));

        if ((projectedEndpoint - projectedVertex).GetMagnitude() > maxPhotonPropagation)
        {
            if (std::find(hierarchy.begin(), hierarchy.end(), pThisMCParticle) != hierarchy.end())
            {
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &projectedEndpoint, "end", BLACK, 2);
                PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &projectedVertex, "start", BLACK, 2);
                std::cout << "Photon Distance: " << (projectedEndpoint - projectedVertex).GetMagnitude() << std::endl;
            }
            return false;
        }
    }

    if (pThisMCParticle == pHitMCParticle)
        return true;

    for (const MCParticle *const pDaughterMCParticle : pThisMCParticle->GetDaughterList())
    {
        if (PrintPhoton(pOriginalPrimary, pDaughterMCParticle, pHitMCParticle, maxPhotonPropagation, hitType))
            return true;
    }

    return false;
}


} // namespace lar_content    
