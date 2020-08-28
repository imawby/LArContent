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
    
    //std::cout << "MCParticle list size: " << pMCParticleList->size() << std::endl;
    //std::cout << "pPfoList size: " << pPfoList->size() << std::endl;
    
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
        LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(*pPfoList, validationInfo.GetAllMCParticleToHitsMap(), pfoToHitsMap, m_primaryParameters.m_foldBackHierarchy);

        validationInfo.SetPfoToHitsMap(pfoToHitsMap);
    }

    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(validationInfo.GetPfoToHitsMap(), {validationInfo.GetAllMCParticleToHitsMap()}, pfoToMCHitSharingMap, mcToPfoHitSharingMap);
    validationInfo.SetMCToPfoHitSharingMap(mcToPfoHitSharingMap);

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMatching(validationInfo, interpretedMCToPfoHitSharingMap);
    validationInfo.SetInterpretedMCToPfoHitSharingMap(interpretedMCToPfoHitSharingMap);
    
    std::cout << "Map size: " << interpretedMCToPfoHitSharingMap.size() << std::endl;
    
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
    
    
    //PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &cosmicRays, "cosmicRays", BLUE);
    //PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &deltaRays, "deltaRays", RED);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayEventValidationAlgorithm::ProcessOutput(const ValidationInfo &/*validationInfo*/, const bool /*useInterpretedMatching*/, const bool /*printToScreen*/, const bool /*fillTree*/) const
{
    
    if (printToScreen && useInterpretedMatching) std::cout << "---INTERPRETED-MATCHING-OUTPUT------------------------------------------------------------------" << std::endl;
    else if (printToScreen) std::cout << "---RAW-MATCHING-OUTPUT--------------------------------------------------------------------------" << std::endl;

    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap(useInterpretedMatching ?
        validationInfo.GetInterpretedMCToPfoHitSharingMap() : validationInfo.GetMCToPfoHitSharingMap());


    typedef std::unordered_map<const pandora::MCParticle*, const MCParticleList> CRToChildDRMap;

    CRToChildDRMap crToChildDRMap;
    
    for (const MCParticleCaloHitListPair &entry : validationInfo.GetAllMCParticleToHitsMap())
    {
        if (!LArMCParticleHelper::IsDeltaRay(entry.first))
            continue;

        const MCParticle *const pParentMuon(LArMCParticleHelper::GetPrimaryMCParticle(entry.first));

        crToChildDRMap[pParentMuon].push_back(entry.first);
    }


    // sort by the highest number of hits CR first 
    MCParticleVector mcCRVector;

    LArMonitoringHelper::GetOrderedMCParticleVector({validationInfo.GetTargetMCParticleToHitsMap()}, mcPrimaryVector);

    // Parent Cosmic Ray 
    int nChildDRs(0), nMatchedChildDRs(0), nMatches_CR(0), isCorrect_CR(0);
    int nMCHitsTotal_CR, nMCHitsU_CR, nMCHitsV_CR, nMCHitsW_CR;
    float mcE_CR, mcPX_CR, mcPY_CR, mcPZ_CR;
    float mcVertexX_CR, mcVertexY_CR, mcVertexZ_CR, mcEndX_CR, mcEndY_CR, mcEndZ_CR;
    float bestMatchNHitsTotal_CR, bestMatchNHitsU_CR, bestMatchNHitsV_CR, bestMatchNHitsW_CR;
    float bestMatchNSharedHitsTotal_CR, bestMatchNSharedHitsU_CR, bestMatchNSharedHitsV_CR, bestMatchNSharedHitsW_CR;
    float bestMatchPfoX0_CR;

    // Child Delta Rays
    IntVector nMatches_DR(0), isCorrect_DR;
    int nMCHitsTotal_DR, nMCHitsU_DR, nMCHitsV_DR, nMCHitsW_DR;
    float mcE_DR, mcPX_DR, mcPY_DR, mcPZ_DR;
    float mcVertexX_DR, mcVertexY_DR, mcVertexZ_DR, mcEndX_DR, mcEndY_DR, mcEndZ_DR;
    float bestMatchNHitsTotal_DR, bestMatchNHitsU_DR, bestMatchNHitsV_DR, bestMatchNHitsW_DR;
    float bestMatchNSharedHitsTotal_DR, bestMatchNSharedHitsU_DR, bestMatchNSharedHitsV_DR, bestMatchNSharedHitsW_DR;
    float bestMatchPfoX0_DR;

    

//
   

}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayEventValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return EventValidationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content    
