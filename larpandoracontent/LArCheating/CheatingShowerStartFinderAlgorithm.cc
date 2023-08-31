/**
 *  @file   larpandoracontent/LArCheating/CheatingShowerStartFinderAlgorithm.cc
 *
 *  @brief  Implementation of the cheating shower start finder
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingShowerStartFinderAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "PandoraMonitoringApi.h"

using namespace pandora;

namespace lar_content
{

CheatingShowerStartFinderAlgorithm::CheatingShowerStartFinderAlgorithm() :
    m_hitDistanceThreshold(1.f),
    m_branchHitThreshold(5),
    m_branchAngleThreshold(20.f),
    m_branchDistanceThreshold(5.f)
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingShowerStartFinderAlgorithm::Run()
{
    // Get MCParticle list
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    // Get CaloHit list
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    // Get the primary electron/photon
    const MCParticle *pMCPrimary(nullptr);

    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsPrimary(pMCParticle))
            pMCPrimary = pMCParticle;
    }

    if (!pMCPrimary)
    {
        std::cout << "No primary MC found!" << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    // Fill CaloHit <-> MCParticle maps
    LArMCParticleHelper::MCRelationMap mcToSelfMap;
    LArMCParticleHelper::CaloHitToMCMap caloHitToMCMap;
    LArMCParticleHelper::MCContributionMap mcToCaloHitListMap; 
    LArMCParticleHelper::GetMCToSelfMap(pMCParticleList, mcToSelfMap);
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcToSelfMap, caloHitToMCMap, mcToCaloHitListMap);

    // Get the true hierarchy of the MC primary
    MCParticleVector particleHierarchy;
    this->FillHierarchyVector(pMCPrimary, particleHierarchy);
    MCParticleVector filteredHierarchyU, filteredHierarchyV, filteredHierarchyW;
    this->PrepareHierarchyVector(pMCPrimary, mcToCaloHitListMap, particleHierarchy, 
        filteredHierarchyU, filteredHierarchyV, filteredHierarchyW);

    const CartesianVector &primaryVertex(pMCPrimary->GetVertex());
    const CartesianVector &primaryDirection(pMCPrimary->GetMomentum().GetUnitVector());
    const CartesianVector primaryEndpoint(primaryVertex + (primaryDirection * 100.0));

    for (HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        MCParticleVector &filteredHierarchyVector(hitType == TPC_VIEW_U ? filteredHierarchyU : hitType == TPC_VIEW_V ? 
            filteredHierarchyV : filteredHierarchyW);

        if (filteredHierarchyVector.empty())
            continue;

        // there's a bug in the simulation...
        CartesianVector viewPrimaryVertex(LArGeometryHelper::ProjectPosition(this->GetPandora(), primaryVertex, hitType));

        float bestXSeparation(std::numeric_limits<float>::max());
        float bestZSeparation(std::numeric_limits<float>::max());
        CartesianVector closestPoint(0.f, 0.f, 0.f);

        for (const CaloHit *const pCaloHit : *pCaloHitList)
        {
            if (pCaloHit->GetHitType() != hitType)
                continue;

            float zSeparation(std::fabs(pCaloHit->GetPositionVector().GetZ() - viewPrimaryVertex.GetZ()));
            float xSeparation(std::fabs(pCaloHit->GetPositionVector().GetX() - viewPrimaryVertex.GetX()));

            if (zSeparation < bestZSeparation)
            {
                bestZSeparation = zSeparation;
                bestXSeparation = xSeparation;
                closestPoint = pCaloHit->GetPositionVector();
            }
        }

        std::cout << "zSeparation: " << bestZSeparation << std::endl;
        std::cout << "xSeparation: " << bestXSeparation << std::endl;

        viewPrimaryVertex = closestPoint;

        const CartesianVector viewPrimaryEndpoint(LArGeometryHelper::ProjectPosition(this->GetPandora(), primaryEndpoint, hitType));
        const MCParticle *const pMCParticle(filteredHierarchyVector.front());
        const CartesianVector viewVertex(LArGeometryHelper::ProjectPosition(this->GetPandora(), pMCParticle->GetVertex(), hitType));

        const CaloHitList &caloHitList(mcToCaloHitListMap.at(pMCParticle));
        CaloHitList viewCaloHitList;
        //float closestDistanceSq(std::numeric_limits<float>::max());
        CartesianVector viewShowerStart(viewVertex);

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            if (pCaloHit->GetHitType() == hitType)
            {
                viewCaloHitList.push_back(pCaloHit);
                //const float separationSq((pCaloHit->GetPositionVector() - viewVertex).GetMagnitudeSquared());

                //if (separationSq < closestDistanceSq)
                //{
                //    closestDistanceSq = separationSq;
                //    viewShowerStart = pCaloHit->GetPositionVector();
                //}
            }
        }

        const CartesianVector viewPrimaryDirection(LArGeometryHelper::ProjectDirection(this->GetPandora(), primaryDirection, hitType));
        const float longitudinalProjection(viewPrimaryDirection.GetDotProduct(viewShowerStart - viewPrimaryVertex));
        viewShowerStart = (viewPrimaryVertex + (viewPrimaryDirection * longitudinalProjection));

        //////////////////////////////////////////////////////////////
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &viewPrimaryVertex, &viewPrimaryEndpoint, "Primary Direction", RED, 2, 2);
        PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &viewCaloHitList, "MCHits", BLACK);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &viewShowerStart, "ShowerStart", BLACK, 2);
        //////////////////////////////////////////////////////////////
    }

    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingShowerStartFinderAlgorithm::FillHierarchyVector(const MCParticle *const pParentMCParticle, MCParticleVector &particleHierarchy) const
{
    for (const MCParticle *const pMCParticle : pParentMCParticle->GetDaughterList())
    {
        if (std::find(particleHierarchy.begin(), particleHierarchy.end(), pMCParticle) != particleHierarchy.end())
            continue;

        particleHierarchy.push_back(pMCParticle);

        this->FillHierarchyVector(pMCParticle, particleHierarchy);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingShowerStartFinderAlgorithm::PrepareHierarchyVector(const MCParticle *const pMCPrimary, 
    const LArMCParticleHelper::MCContributionMap &mcToCaloHitListMap, const MCParticleVector &particleHierarchy, 
    MCParticleVector &filteredHierarchyU, MCParticleVector &filteredHierarchyV, MCParticleVector &filteredHierarchyW) const
{
    const CartesianVector &primaryVertex(pMCPrimary->GetVertex());
    const CartesianVector primaryVertexU(LArGeometryHelper::ProjectPosition(this->GetPandora(), primaryVertex, TPC_VIEW_U));
    const CartesianVector primaryVertexV(LArGeometryHelper::ProjectPosition(this->GetPandora(), primaryVertex, TPC_VIEW_V));
    const CartesianVector primaryVertexW(LArGeometryHelper::ProjectPosition(this->GetPandora(), primaryVertex, TPC_VIEW_W));
    //const CartesianVector &primaryEndpoint(pMCPrimary->GetEndpoint());
    //const CartesianVector primaryEndpointU(LArGeometryHelper::ProjectPosition(this->GetPandora(), primaryEndpoint, TPC_VIEW_U));
    //const CartesianVector primaryEndpointV(LArGeometryHelper::ProjectPosition(this->GetPandora(), primaryEndpoint, TPC_VIEW_V));
    //const CartesianVector primaryEndpointW(LArGeometryHelper::ProjectPosition(this->GetPandora(), primaryEndpoint, TPC_VIEW_W));
    const CartesianVector &primaryDirection(pMCPrimary->GetMomentum().GetUnitVector());
    //const CartesianVector primaryDirectionU((primaryEndpointU - primaryVertexU).GetUnitVector());
    //const CartesianVector primaryDirectionV((primaryEndpointV - primaryVertexV).GetUnitVector());
    //const CartesianVector primaryDirectionW((primaryEndpointW - primaryVertexW).GetUnitVector());
    const CartesianVector primaryDirectionU(LArGeometryHelper::ProjectDirection(this->GetPandora(), primaryDirection, TPC_VIEW_U));
    const CartesianVector primaryDirectionV(LArGeometryHelper::ProjectDirection(this->GetPandora(), primaryDirection, TPC_VIEW_V));
    const CartesianVector primaryDirectionW(LArGeometryHelper::ProjectDirection(this->GetPandora(), primaryDirection, TPC_VIEW_W));

    for (const MCParticle *const pMCParticle : particleHierarchy)
    {
        if (mcToCaloHitListMap.find(pMCParticle) == mcToCaloHitListMap.end())
            continue;

        // Collect hits
        const CaloHitList &caloHitList(mcToCaloHitListMap.at(pMCParticle));
        unsigned int caloHitCountU(0), caloHitCountV(0), caloHitCountW(0);

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            const HitType hitType(pCaloHit->GetHitType());
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            const CartesianVector &primaryVertex2D(hitType == TPC_VIEW_U ? primaryVertexU : hitType == TPC_VIEW_V ? 
                primaryVertexV : primaryVertexW);
            const CartesianVector &primaryDirection2D(hitType == TPC_VIEW_U ? primaryDirectionU : hitType == TPC_VIEW_V ? 
                primaryDirectionV : primaryDirectionW);
            const float separationSquared(primaryDirection2D.GetCrossProduct(hitPosition - primaryVertex2D).GetMagnitudeSquared());

            if ((pMCPrimary->GetParticleId() == 22) && (separationSquared < (m_hitDistanceThreshold * m_hitDistanceThreshold)))
                continue;

            unsigned int &caloHitCount(hitType == TPC_VIEW_U ? caloHitCountU : hitType == TPC_VIEW_V ? caloHitCountV : caloHitCountW);
            ++caloHitCount;
        }

        // Calculate deviation threshold
        const CartesianVector &direction(pMCParticle->GetMomentum().GetUnitVector());
        const float openingAngle(direction.GetOpeningAngle(primaryDirection) * 180.0 / 3.14);

        // Apply degree threshold
        if (openingAngle < m_branchAngleThreshold)
           continue;

        // Calculate transverse separation distance
        const CartesianVector &vertex(pMCParticle->GetVertex());
        const float transverseDistance(direction.GetCrossProduct(vertex - primaryVertex).GetMagnitude());

        // Apply distance threshold
        if (transverseDistance < m_branchDistanceThreshold)
            continue;

        if (caloHitCountU >= m_branchHitThreshold)
            filteredHierarchyU.push_back(pMCParticle);

        if (caloHitCountV >= m_branchHitThreshold)
            filteredHierarchyV.push_back(pMCParticle);

        if (caloHitCountW >= m_branchHitThreshold)
            filteredHierarchyW.push_back(pMCParticle);
    }

    for (HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        MCParticleVector &filteredHierarchyVector(hitType == TPC_VIEW_U ? filteredHierarchyU : hitType == TPC_VIEW_V ? 
            filteredHierarchyV : filteredHierarchyW);

        std::sort(filteredHierarchyVector.begin(), filteredHierarchyVector.end(), 
            [primaryVertex](const MCParticle *const &pMCParticleA, const MCParticle *const pMCParticleB) -> bool
            { 
                const float separationSquaredA((primaryVertex - pMCParticleA->GetVertex()).GetMagnitudeSquared());
                const float separationSquaredB((primaryVertex - pMCParticleB->GetVertex()).GetMagnitudeSquared());

                return separationSquaredA < separationSquaredB;
            }
        );
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingShowerStartFinderAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "HitDistanceThreshold", m_hitDistanceThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "BranchHitThreshold", m_branchHitThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "BranchAngleThreshold", m_branchAngleThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "BranchDistanceThreshold", m_branchDistanceThreshold));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
