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
    
    std::cout << "MCParticle list size: " << pMCParticleList->size() << std::endl;
    std::cout << "pPfoList size: " << pPfoList->size() << std::endl;

    /*
    CartesianPointVector deltaRayHits, cosmicRayHits;
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        const MCParticle *pMCParticle(nullptr);
        
        try
        {
            pMCParticle = MCParticleHelper::GetMainMCParticle(pCaloHit);
        }
        catch (StatusCodeException&)
        {
        }

        if (pMCParticle)
        {
            if ( LArMCParticleHelper::IsCosmicRay(pMCParticle) && LArMCParticleHelper::IsDeltaRay(pMCParticle))
            {
                std::cout << " hit belongs to both a cosmic ray and a delta ray!" << std::endl;
                continue;
            }

            if (LArMCParticleHelper::IsCosmicRay(pMCParticle))
                cosmicRayHits.push_back(pCaloHit->GetPositionVector());

            if (LArMCParticleHelper::IsDeltaRay(pMCParticle))
                deltaRayHits.push_back(pCaloHit->GetPositionVector());

        }
    }

    for (const CartesianVector &point : deltaRayHits)
    {
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &point, "DR", RED, 2);
    }

    for (const CartesianVector &point : cosmicRayHits)
    {
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &point, "CR", BLUE, 2);
    }

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    
    
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        std::cout << "calo hit map size: " << pCaloHit->GetMCParticleWeightMap().size() << std::endl;
    }

    PandoraMonitoringApi::Create(this->GetPandora());
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);
    
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

    LArMCParticleHelper::MCParticleToPfoHitSharingMap interpretedMCToPfoHitSharingMap;
    this->InterpretMatching(validationInfo, interpretedMCToPfoHitSharingMap);
    validationInfo.SetInterpretedMCToPfoHitSharingMap(interpretedMCToPfoHitSharingMap);
    
    std::cout << "Map size: " << interpretedMCToPfoHitSharingMap.size() << std::endl;
    
    PfoList cosmicRays, deltaRays;
    unsigned int count_DR(0), hasMatch_DR(0);
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
            count_DR++;

            if (!entry.second.empty())
            {
                hasMatch_DR++;
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

        const ParticleFlowObject *pParentMuon(LArPfoHelper::GetParentPfo(pDeltaRay));

        
        CartesianPointVector position;
        LArPfoHelper::GetCoordinateVector(pDeltaRay, TPC_3D, position);

        for (const CartesianVector &point : position)
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &point, "DR", RED, 2);
    }    
    

    std::cout << "count_DR: " << count_DR << std::endl;
    std::cout << "hasMatch_DR: " << hasMatch_DR << std::endl;
    
    //PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &cosmicRays, "cosmicRays", BLUE);
    //PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &deltaRays, "deltaRays", RED);
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
