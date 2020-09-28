/**
 *  @file   larpandoracontent/LArMonitoring/DeltaRayEventValidationAlgorithm.cc
 *
 *  @brief  Implementation of the delta ray event validation algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/DeltaRayEventValidationAlgorithm.h"

#include <sstream>

using namespace pandora;

namespace lar_content
{

DeltaRayEventValidationAlgorithm::DeltaRayEventValidationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DeltaRayEventValidationAlgorithm::~DeltaRayEventValidationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayEventValidationAlgorithm::FillValidationInfo(const MCParticleList *const pMCParticleList, const CaloHitList *const pCaloHitList,
    const PfoList *const pPfoList, ValidationInfo &validationInfo) const
{
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);
    
    if (pMCParticleList && pCaloHitList)
    {
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, m_primaryParameters, LArMCParticleHelper::IsCosmicRay, targetMCParticleToHitsMap);

        LArMCParticleHelper::PrimaryParameters parameters(m_primaryParameters);
        parameters.m_minPrimaryGoodHits = 0;
        parameters.m_minHitsForGoodView = 0;
        parameters.m_minHitSharingFraction = 0.f;
        LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, allMCParticleToHitsMap);

        validationInfo.SetTargetMCParticleToHitsMap(targetMCParticleToHitsMap);
        validationInfo.SetAllMCParticleToHitsMap(allMCParticleToHitsMap);
    }

    if (pPfoList)
    {
        LArMCParticleHelper::PfoContributionMap pfoToHitsMap;
        //LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(*pPfoList, validationInfo.GetAllMCParticleToHitsMap(), pfoToHitsMap, m_primaryParameters.m_foldBackHierarchy);
        LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(*pPfoList, validationInfo.GetTargetMCParticleToHitsMap(), pfoToHitsMap, m_primaryParameters.m_foldBackHierarchy);

        validationInfo.SetPfoToHitsMap(pfoToHitsMap);
    }

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    //LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(validationInfo.GetPfoToHitsMap(), {validationInfo.GetAllMCParticleToHitsMap()}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(validationInfo.GetPfoToHitsMap(), {validationInfo.GetTargetMCParticleToHitsMap()}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);
    
    validationInfo.SetMCToPfoHitSharingMap(mcToPfoHitSharingMap);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMatching(validationInfo, interpretedMCToPfoHitSharingMap);
    validationInfo.SetInterpretedMCToPfoHitSharingMap(interpretedMCToPfoHitSharingMap);

    std::cout << "MC Hit Map size: " << validationInfo.GetTargetMCParticleToHitsMap().size() << std::endl;
    std::cout << "Pfo Hit Map size: " << validationInfo.GetPfoToHitsMap().size() << std::endl;
    std::cout << "Matching Map size: " << mcToPfoHitSharingMap.size() << std::endl;
    std::cout << "Interpretted Matching Map size: " << interpretedMCToPfoHitSharingMap.size() << std::endl;

    /*
    for (auto &entry : validationInfo.GetTargetMCParticleToHitsMap())
    {
        if (mcToPfoHitSharingMap.find(entry.first) == mcToPfoHitSharingMap.end())
            continue;

        std::cout << "HITS: " << entry.second.size() << std::endl;
        
    }
    */
    
    PfoList cosmicRays, deltaRays;
    for (const auto &entry : interpretedMCToPfoHitSharingMap)
    {
        if (LArMCParticleHelper::IsCosmicRay(entry.first))
        {
            if (!entry.second.empty())
            {
                cosmicRays.push_back(entry.second.begin()->first);
            }
        }

        if (LArMCParticleHelper::IsDeltaRay(entry.first))
        {
            if (!entry.second.empty())
            {
                deltaRays.push_back(entry.second.begin()->first);
            }
        }
    }

    for (const ParticleFlowObject *const pCosmicRay : cosmicRays)
    {
        CartesianPointVector position;
        LArPfoHelper::GetCoordinateVector(pCosmicRay, TPC_3D, position);

        for (const CartesianVector &point : position)
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &point, "CR", BLUE, 2);
    }

    for (const ParticleFlowObject *const pDeltaRay : deltaRays)
    {
        CartesianPointVector position;
        LArPfoHelper::GetCoordinateVector(pDeltaRay, TPC_3D, position);   

        for (const CartesianVector &point : position)
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &point, "DR", RED, 2);
    }    
    
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayEventValidationAlgorithm::ProcessOutput(const ValidationInfo &validationInfo, const bool useInterpretedMatching, const bool printToScreen, const bool fillTree) const
{
    if (printToScreen && useInterpretedMatching) std::cout << "---INTERPRETED-MATCHING-OUTPUT------------------------------------------------------------------" << std::endl;
    else if (printToScreen) std::cout << "---RAW-MATCHING-OUTPUT--------------------------------------------------------------------------" << std::endl;

    
    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap(useInterpretedMatching ?
        validationInfo.GetInterpretedMCToPfoHitSharingMap() : validationInfo.GetMCToPfoHitSharingMap());
    
    LArMCParticleHelper::CRToChildDRMap crToChildDRMap;
    const LArMCParticleHelper::MCContributionMap &targetMCParticleToHitsMap(validationInfo.GetTargetMCParticleToHitsMap());
    for (const LArMCParticleHelper::MCParticleCaloHitListPair &entry : targetMCParticleToHitsMap)
    {
        if (!LArMCParticleHelper::IsDeltaRay(entry.first))
            continue;

        const MCParticle *const pParentMuon(LArMCParticleHelper::GetPrimaryMCParticle(entry.first));

        if (targetMCParticleToHitsMap.find(pParentMuon) == targetMCParticleToHitsMap.end())
            continue;

        crToChildDRMap[pParentMuon].push_back(entry.first);
    }
    
    // sort by the highest number of hits CR first 
    MCParticleVector mcCRVector;
    for (auto &entry : crToChildDRMap)
        mcCRVector.push_back(entry.first);
    
    LArMonitoringHelper::GetOrderedMCParticleVector(targetMCParticleToHitsMap, mcCRVector);
  
    // Parent Cosmic Ray 
    int nReconstructableChildDRs, nMatches_CR(0), isCorrect_CR, nCorrectChildDRs(0);
    float mcE_CR, mcPX_CR, mcPY_CR, mcPZ_CR;    
    int nMCHitsTotal_CR, nMCHitsU_CR, nMCHitsV_CR, nMCHitsW_CR;
    float mcVertexX_CR, mcVertexY_CR, mcVertexZ_CR, mcEndX_CR, mcEndY_CR, mcEndZ_CR;
    float bestMatchNHitsTotal_CR(0), bestMatchNHitsU_CR(0), bestMatchNHitsV_CR(0), bestMatchNHitsW_CR(0);
    float bestMatchNSharedHitsTotal_CR(0), bestMatchNSharedHitsU_CR(0), bestMatchNSharedHitsV_CR(0), bestMatchNSharedHitsW_CR(0);
     //float bestMatchPfoX0_CR;
    
    // Child Delta Rays
    IntVector nMatches_DR, isCorrect_DR;
    FloatVector mcE_DR, mcPX_DR, mcPY_DR, mcPZ_DR;    
    IntVector nMCHitsTotal_DR, nMCHitsU_DR, nMCHitsV_DR, nMCHitsW_DR;
    FloatVector mcVertexX_DR, mcVertexY_DR, mcVertexZ_DR, mcEndX_DR, mcEndY_DR, mcEndZ_DR;
    FloatVector bestMatchNHitsTotal_DR, bestMatchNHitsU_DR, bestMatchNHitsV_DR, bestMatchNHitsW_DR;
    FloatVector bestMatchNSharedHitsTotal_DR, bestMatchNSharedHitsU_DR, bestMatchNSharedHitsV_DR, bestMatchNSharedHitsW_DR;
    //float bestMatchPfoX0_DR;

    std::cout << "HERE" << std::endl;
    
    std::stringstream stringStream;
    for (const MCParticle *const pCosmicRay : mcCRVector)
    {
        const CaloHitList &cosmicRayHitList(validationInfo.GetAllMCParticleToHitsMap().at(pCosmicRay));

        std::cout << "AAAAAAAAAAA" << std::endl;
        nReconstructableChildDRs = crToChildDRMap.at(pCosmicRay).size();
        std::cout << "BBBBBBBBBBB" << std::endl;
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

        // Add angle?
        stringStream << "(Parent CR) "
            << "Energy " << pCosmicRay->GetEnergy()
            << ", Dist. " << (pCosmicRay->GetEndpoint() - pCosmicRay->GetVertex()).GetMagnitude()
            << ", nMCHits " << cosmicRayHitList.size()
            << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, cosmicRayHitList)
            << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, cosmicRayHitList)
            << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, cosmicRayHitList) << ")" << std::endl;

        std::cout << "CCCCCCCCCCCC" << std::endl;
        for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcToPfoHitSharingMap.at(pCosmicRay))
        {
            
            std::cout << "111111111111111111111111111111111" << std::endl;
            
            const CaloHitList &sharedHitList(pfoToSharedHits.second);
            const CaloHitList &pfoHitList(validationInfo.GetPfoToHitsMap().at(pfoToSharedHits.first));

            const bool isGoodMatch_CR(this->IsGoodMatch(cosmicRayHitList, pfoHitList, sharedHitList)); 

            if (0 == nMatches_CR++)
            {
                bestMatchNHitsTotal_CR = pfoHitList.size();
                bestMatchNHitsU_CR = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList);
                bestMatchNHitsV_CR = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList);
                bestMatchNHitsW_CR = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList);
                bestMatchNSharedHitsTotal_CR = sharedHitList.size();
                bestMatchNSharedHitsU_CR = LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList);
                bestMatchNSharedHitsV_CR = LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList);
                bestMatchNSharedHitsW_CR = LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList);
            }

            std::cout << "DDDDDDDDDDDDDDDDD" << std::endl;

            stringStream << "-" << (!isGoodMatch_CR ? "(Below threshold) " : "")
                << "nMatchedHits " << sharedHitList.size()
                << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList)
                << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList)
                << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList) << ")"
                << ", nPfoHits " << pfoHitList.size()
                << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList)
                << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList)
                << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList) << ")" << std::endl;

             if (isGoodMatch_CR)
                 ++nMatches_CR;

             std::cout << "EEEEEEEEEEEEEEEE" << std::endl;
        }

        if (mcToPfoHitSharingMap.at(pCosmicRay).empty())
            stringStream << "-No matched Pfo" << std::endl;

        isCorrect_CR = (nMatches_CR == 1);

        stringStream << (isCorrect_CR ? "Correct CR" : "Incorrect CR") << std::endl;

        // Look at child delta ray particles
        std::cout << "AAAAAAAAAAA" << std::endl;
        MCParticleVector childDeltaRayMCParticles(crToChildDRMap.at(pCosmicRay));
        std::cout << "BBBBBBBBBBB" << std::endl;
        LArMonitoringHelper::GetOrderedMCParticleVector(targetMCParticleToHitsMap, childDeltaRayMCParticles);
        for (const MCParticle *const pDeltaRay : childDeltaRayMCParticles)
        {
            const CaloHitList &deltaRayHitList(validationInfo.GetAllMCParticleToHitsMap().at(pDeltaRay));

            mcE_DR.push_back(pCosmicRay->GetEnergy());
            mcPX_DR.push_back(pCosmicRay->GetMomentum().GetX());
            mcPY_DR.push_back(pCosmicRay->GetMomentum().GetY());
            mcPZ_DR.push_back(pCosmicRay->GetMomentum().GetZ());
            mcVertexX_DR.push_back(pCosmicRay->GetVertex().GetX());
            mcVertexY_DR.push_back(pCosmicRay->GetVertex().GetY());
            mcVertexZ_DR.push_back(pCosmicRay->GetVertex().GetZ());        
            mcEndX_DR.push_back(pCosmicRay->GetEndpoint().GetX());
            mcEndY_DR.push_back(pCosmicRay->GetEndpoint().GetY());
            mcEndZ_DR.push_back(pCosmicRay->GetEndpoint().GetZ());   
            nMCHitsTotal_DR.push_back(cosmicRayHitList.size());
            nMCHitsU_DR.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, cosmicRayHitList));
            nMCHitsV_DR.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, cosmicRayHitList));   
            nMCHitsW_DR.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, cosmicRayHitList));

            // Add angle?
            stringStream << "(Child DR) "
                << "Energy " << pDeltaRay->GetEnergy()
                << ", Dist. " << (pDeltaRay->GetEndpoint() - pDeltaRay->GetVertex()).GetMagnitude()
                << ", nMCHits " << deltaRayHitList.size()
                << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, deltaRayHitList)
                << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, deltaRayHitList)
                 << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, deltaRayHitList) << ")" << std::endl;

            int nMatches(0);
            for (const LArMCParticleHelper::PfoCaloHitListPair &pfoToSharedHits : mcToPfoHitSharingMap.at(pDeltaRay))
            {
                const CaloHitList &sharedHitList(pfoToSharedHits.second);
                const CaloHitList &pfoHitList(validationInfo.GetPfoToHitsMap().at(pfoToSharedHits.first));

                const bool isGoodMatch_DR(this->IsGoodMatch(deltaRayHitList, pfoHitList, sharedHitList));

                if (0 == nMatches++)
                {
                    bestMatchNHitsTotal_DR.push_back(pfoHitList.size());
                    bestMatchNHitsU_DR.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList));
                    bestMatchNHitsV_DR.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList));
                    bestMatchNHitsW_DR.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList));
                    bestMatchNSharedHitsTotal_DR.push_back(sharedHitList.size());
                    bestMatchNSharedHitsU_DR.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList));
                    bestMatchNSharedHitsV_DR.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList));
                    bestMatchNSharedHitsW_DR.push_back(LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList));
                }

                stringStream << "-" << (!isGoodMatch_DR ? "(Below threshold) " : "")
                    << "nMatchedHits " << sharedHitList.size()
                    << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, sharedHitList)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, sharedHitList)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, sharedHitList) << ")"
                    << ", nPfoHits " << pfoHitList.size()
                    << " (" << LArMonitoringHelper::CountHitsByType(TPC_VIEW_U, pfoHitList)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_V, pfoHitList)
                    << ", " << LArMonitoringHelper::CountHitsByType(TPC_VIEW_W, pfoHitList) << ")" << std::endl;

                if (isGoodMatch_DR)
                    ++nMatches;
            }

            stringStream << (nMatches == 1 ? "Correct DR" : "Incorrect DR") << std::endl;
            
            if (nMatches == 1)
            {
                isCorrect_DR.push_back(1);
                ++nCorrectChildDRs;
            }

            if (mcToPfoHitSharingMap.at(pDeltaRay).empty())
            {
                stringStream << "-No matched Pfo" << std::endl;
                bestMatchNHitsTotal_DR.push_back(0); bestMatchNHitsU_DR.push_back(0); bestMatchNHitsV_DR.push_back(0); bestMatchNHitsW_DR.push_back(0);
                bestMatchNSharedHitsTotal_DR.push_back(0); bestMatchNSharedHitsU_DR.push_back(0); bestMatchNSharedHitsV_DR.push_back(0); bestMatchNSharedHitsW_DR.push_back(0);
                nMatches_DR.push_back(0), isCorrect_DR.push_back(0);
            }
        }

	    if (fillTree)
        {
            // Write general event variables
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "fileIdentifier", m_fileIdentifier));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "eventNumber", m_eventNumber - 1));

            // Write CR variables
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nReconstructableChildDRs", nReconstructableChildDRs));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMatches_CR", nMatches_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrect_CR", isCorrect_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nCorrectChildDRs", nCorrectChildDRs));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsTotal_CR", nMCHitsTotal_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsU_CR", nMCHitsU_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsV_CR", nMCHitsV_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsW_CR", nMCHitsW_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcE_CR", mcE_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPX_CR", mcPX_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPY_CR", mcPY_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPZ_CR", mcPZ_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexX_CR", mcVertexX_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexY_CR", mcVertexY_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexZ_CR", mcVertexZ_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndX_CR", mcEndX_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndY_CR", mcEndY_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndZ_CR", mcEndZ_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNHitsTotal_CR", bestMatchNHitsTotal_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNHitsU_CR", bestMatchNHitsU_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNHitsV_CR", bestMatchNHitsV_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNHitsW_CR", bestMatchNHitsW_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNSharedHitsTotal_CR", bestMatchNSharedHitsTotal_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNSharedHitsU_CR", bestMatchNSharedHitsU_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNSharedHitsV_CR", bestMatchNSharedHitsV_CR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNSharedHitsW_CR", bestMatchNSharedHitsW_CR));

            // Write DR variables
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMatches_DR", &nMatches_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "isCorrect_DR", &isCorrect_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsTotal_DR", &nMCHitsTotal_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsU_DR", &nMCHitsU_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsV_DR", &nMCHitsV_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "nMCHitsW_DR", &nMCHitsW_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcE_DR", &mcE_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPX_DR", &mcPX_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPY_DR", &mcPY_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcPZ_DR", &mcPZ_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexX_DR", &mcVertexX_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexY_DR", &mcVertexY_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcVertexZ_DR", &mcVertexZ_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndX_DR", &mcEndX_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndY_DR", &mcEndY_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "mcEndZ_DR", &mcEndZ_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNHitsTotal_DR", &bestMatchNHitsTotal_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNHitsU_DR", &bestMatchNHitsU_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNHitsV_DR", &bestMatchNHitsV_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNHitsW_DR", &bestMatchNHitsW_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNSharedHitsTotal_DR", &bestMatchNSharedHitsTotal_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNSharedHitsU_DR", &bestMatchNSharedHitsU_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNSharedHitsV_DR", &bestMatchNSharedHitsV_DR));
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName.c_str(), "bestMatchNSharedHitsW_DR", &bestMatchNSharedHitsW_DR));
        }

        // Reset variables
        nMatches_CR = 0; nCorrectChildDRs = 0;
        bestMatchNHitsTotal_CR = 0; bestMatchNHitsU_CR = 0; bestMatchNHitsV_CR = 0; bestMatchNHitsW_CR = 0;
        bestMatchNSharedHitsTotal_CR = 0; bestMatchNSharedHitsU_CR = 0; bestMatchNSharedHitsV_CR = 0; bestMatchNSharedHitsW_CR = 0;
        
        nMatches_DR.clear(); isCorrect_DR.clear();
        mcE_DR.clear(); mcPX_DR.clear(); mcPY_DR.clear(); mcPZ_DR.clear();
        nMCHitsTotal_DR.clear(); nMCHitsU_DR.clear(); nMCHitsV_DR.clear(); nMCHitsW_DR.clear();

        mcVertexX_DR.clear(); mcVertexY_DR.clear(); mcVertexZ_DR.clear(); mcEndX_DR.clear(); mcEndY_DR.clear(); mcEndZ_DR.clear();
        bestMatchNHitsTotal_DR.clear(); bestMatchNHitsU_DR.clear(); bestMatchNHitsV_DR.clear(); bestMatchNHitsW_DR.clear();
        bestMatchNSharedHitsTotal_DR.clear(); bestMatchNSharedHitsU_DR.clear(); bestMatchNSharedHitsV_DR.clear(); bestMatchNSharedHitsW_DR.clear();

        if (printToScreen)
            std::cout << stringStream.str() << std::endl;

        if (printToScreen) std::cout << "------------------------------------------------------------------------------------------------" << std::endl << std::endl;
    }


        

  


    

}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayEventValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return EventValidationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content    
