/**
 *  @file   larpandoracontent/LArMonitoring/MCDeltaRay.cc
 *
 *  @brief  Implementation of the mc delta ray algorithm.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "larpandoracontent/LArMonitoring/MCDeltaRayAlgorithm.h"

#include "larpandoracontent/LArObjects/LArMCParticle.h"

using namespace pandora;

namespace lar_content
{

MCDeltaRayAlgorithm::MCDeltaRayAlgorithm() :
    m_caloHitListName(),
    m_mcParticleListName()
{    
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCDeltaRayAlgorithm::Run()
{

    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);
    
    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));





    
    std::cout << "MCParticle list size: " << pMCParticleList->size() << std::endl;

    LArMCParticleHelper::PrimaryParameters parameters;
    parameters.m_foldBackHierarchy = false;

    LArMCParticleHelper::MCContributionMap targetMCParticleToHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, parameters, LArMCParticleHelper::IsCosmicRay, targetMCParticleToHitsMap);


   
    for (auto &entry : targetMCParticleToHitsMap)
    {
        const MCParticle *const pMCParticle(entry.first);
        
        if (!LArMCParticleHelper::IsDeltaRay(pMCParticle))
            continue;
        
        const MCParticle *const pParentMuon(LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle));
        auto parentMuonIter(targetMCParticleToHitsMap.find(pParentMuon));

        std::cout << "Parent 'muon' PDG: " << std::fabs(parentMuonIter->first->GetParticleId()) << std::endl;

        if (parentMuonIter == targetMCParticleToHitsMap.end())
            continue;

        CartesianPointVector deltaRayHits, cosmicRayHits;
        CartesianPointVector deltaRayHits_U, deltaRayHits_V, deltaRayHits_W;
        CartesianPointVector cosmicRayHits_U, cosmicRayHits_V, cosmicRayHits_W;
        
        for (const CaloHit *const pCaloHit : entry.second)
        {
            deltaRayHits.push_back(pCaloHit->GetPositionVector());

            if (pCaloHit->GetHitType() == TPC_VIEW_U)
                deltaRayHits_U.push_back(pCaloHit->GetPositionVector());

            if (pCaloHit->GetHitType() == TPC_VIEW_V)
                deltaRayHits_V.push_back(pCaloHit->GetPositionVector());

            if (pCaloHit->GetHitType() == TPC_VIEW_W)
                deltaRayHits_W.push_back(pCaloHit->GetPositionVector());            
        }

        for (const CaloHit *const pCaloHit : parentMuonIter->second)
        {
            cosmicRayHits.push_back(pCaloHit->GetPositionVector());

            if (pCaloHit->GetHitType() == TPC_VIEW_U)
                cosmicRayHits_U.push_back(pCaloHit->GetPositionVector());

            if (pCaloHit->GetHitType() == TPC_VIEW_V)
                cosmicRayHits_V.push_back(pCaloHit->GetPositionVector());

            if (pCaloHit->GetHitType() == TPC_VIEW_W)
                cosmicRayHits_W.push_back(pCaloHit->GetPositionVector());            
        }

        for (const CartesianVector &point : deltaRayHits)
        {
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &point, "DR", RED, 1);
        }

        for (const CartesianVector &point : cosmicRayHits)
        {
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &point, "CR", BLUE, 1);
        }

        CartesianVector muonVertex_U(0.f, 0.f, 0.f), muonEndpoint_U(0.f, 0.f, 0.f);
        LArClusterHelper::GetExtremalCoordinates(cosmicRayHits_U, muonEndpoint_U, muonVertex_U);

        CartesianVector muonVertex_V(0.f, 0.f, 0.f), muonEndpoint_V(0.f, 0.f, 0.f);
        LArClusterHelper::GetExtremalCoordinates(cosmicRayHits_V, muonEndpoint_V, muonVertex_V);

        CartesianVector muonVertex_W(0.f, 0.f, 0.f), muonEndpoint_W(0.f, 0.f, 0.f);
        LArClusterHelper::GetExtremalCoordinates(cosmicRayHits_W, muonEndpoint_W, muonVertex_W);        

        float muonVertexChiSquared(0.f);
        CartesianVector muonVertex3D(0.f, 0.f, 0.f);
        LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, muonVertex_U, muonVertex_V, muonVertex_W, muonVertex3D, muonVertexChiSquared);

        float muonEndpointChiSquared(0.f);
        CartesianVector muonEndpoint3D(0.f, 0.f, 0.f);
        LArGeometryHelper::MergeThreePositions3D(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W, muonEndpoint_U, muonEndpoint_V, muonEndpoint_W, muonEndpoint3D, muonEndpointChiSquared);
     
        
        const CartesianVector &mcParticleVertex(pMCParticle->GetVertex());

        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &muonEndpoint3D, "muonEndpoint", BLACK, 3);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &muonVertex3D, "muonVertex", BLACK, 3);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &mcParticleVertex, "mcParticleVertex", BLACK, 3);

        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }



    

    /*
    CartesianPointVector deltaRayHits, cosmicRayHits;
    for (auto &entry : targetMCParticleToHitsMap)
    {
        
        for (const CaloHit *const pCaloHit : entry.second)
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
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &point, "DR", RED, 1);
    }

    for (const CartesianVector &point : cosmicRayHits)
    {
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &point, "CR", BLUE, 1);
    }

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    
    return STATUS_CODE_SUCCESS;


}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode MCDeltaRayAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
