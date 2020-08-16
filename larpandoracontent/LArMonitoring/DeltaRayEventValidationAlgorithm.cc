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

    std::cout << "MCParticle list size: " << pMCParticleList->size() << std::endl;
    
    if (pMCParticleList && pCaloHitList)
    {
        LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, m_primaryParameters, LArMCParticleHelper::IsCosmicRay, targetMCParticleToHitsMap);
        //LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, m_primaryParameters, LArMCParticleHelper::IsDeltaRay, targetMCParticleToHitsMap);        

        LArMCParticleHelper::PrimaryParameters parameters(m_primaryParameters);
        parameters.m_minPrimaryGoodHits = 0;
        parameters.m_minHitsForGoodView = 0;
        parameters.m_minHitSharingFraction = 0.f;
        LArMCParticleHelper::MCContributionMap allMCParticleToHitsMap;
        LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, allMCParticleToHitsMap);
        //LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsDeltaRay, allMCParticleToHitsMap);        

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

    PfoList cosmicRays, deltaRays;
    for (const auto &entry : mcToPfoHitSharingMap)
    {
        if (LArMCParticleHelper::IsCosmicRay(entry.first))
            cosmicRays.push_back(entry.second.begin()->first);

        if (LArMCParticleHelper::IsDeltaRay(entry.first))
            deltaRays.push_back(entry.second.begin()->first);
    }

    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &cosmicRays, "cosmicRays", BLUE);
    PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &deltaRays, "deltaRays", RED);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DeltaRayEventValidationAlgorithm::ProcessOutput(const ValidationInfo &/*validationInfo*/, const bool /*useInterpretedMatching*/, const bool /*printToScreen*/, const bool /*fillTree*/) const
{
    /*
    if (printToScreen && useInterpretedMatching) std::cout << "---INTERPRETED-MATCHING-OUTPUT------------------------------------------------------------------" << std::endl;
    else if (printToScreen) std::cout << "---RAW-MATCHING-OUTPUT--------------------------------------------------------------------------" << std::endl;

    const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcToPfoHitSharingMap(useInterpretedMatching ?
        validationInfo.GetInterpretedMCToPfoHitSharingMap() : validationInfo.GetMCToPfoHitSharingMap());

    MCParticleVector mcPrimaryVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({validationInfo.GetTargetMCParticleToHitsMap()}, mcPrimaryVector);

    int correctParentCR(0), correctDR(0);
    int nParentCRMatches(0), nDRMatches);
    IntVector mcPrimaryId, mcPrimaryPdg, nMCHitsTotal, nMCHitsU, nMCHitsV, nMCHitsW;
    FloatVector mcPrimaryE, mcPrimaryPX, mcPrimaryPY, mcPrimaryPZ;
    FloatVector mcPrimaryVtxX, mcPrimaryVtxY, mcPrimaryVtxZ, mcPrimaryEndX, mcPrimaryEndY, mcPrimaryEndZ;
    IntVector nPrimaryMatchedPfos, nPrimaryMatchedTBPfos, nPrimaryMatchedCRPfos;
    IntVector bestMatchPfoId, bestMatchPfoPdg, bestMatchPfoIsTB;
    IntVector bestMatchPfoNHitsTotal, bestMatchPfoNHitsU, bestMatchPfoNHitsV, bestMatchPfoNHitsW;
    IntVector bestMatchPfoNSharedHitsTotal, bestMatchPfoNSharedHitsU, bestMatchPfoNSharedHitsV, bestMatchPfoNSharedHitsW;
    FloatVector bestMatchPfoX0;
    */

}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DeltaRayEventValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return EventValidationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content    
