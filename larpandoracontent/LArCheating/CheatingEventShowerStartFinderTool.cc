/**
 *  @file   larpandoracontent/LArShowerRefinement/CheatingEventShowerStartFinderTool.cc
 *
 *  @brief  Implementation of the shower start finder tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"


#include "larpandoracontent/LArCheating/CheatingEventShowerStartFinderTool.h"

using namespace pandora;

namespace lar_content
{

CheatingEventShowerStartFinderTool::CheatingEventShowerStartFinderTool() : 
    m_threshold2DHitCount(100),
    m_hitDistanceThreshold(1.f),
    m_branchHitThreshold(5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingEventShowerStartFinderTool::Run(const MCParticleList *const pMCParticleList, const CaloHitList *const pCaloHitList, 
    CartesianPointVector &showerStartsU, CartesianPointVector &showerStartsV, CartesianPointVector &showerStartsW)
{
    // This tool should run once per-event
    this->ResetMaps();
    this->FillTrueHitMaps(pMCParticleList, pCaloHitList);

    // Get all electrons/photons that are chunky?
    TargetHierarchyMap targetHierarchyMap;
    this->GetTargetShowers(pMCParticleList, targetHierarchyMap);

    for (const auto &entry : targetHierarchyMap)
    {
        for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
        {
            CartesianVector viewShowerStart(0.f, 0.f, 0.f);
            CartesianPointVector &viewShowerStarts(hitType == TPC_VIEW_U ? showerStartsU : hitType == TPC_VIEW_V ? showerStartsV : showerStartsW);

            if (this->GetViewShowerStart(entry.first, entry.second, hitType, viewShowerStart))
                viewShowerStarts.push_back(viewShowerStart);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingEventShowerStartFinderTool::ResetMaps()
{
    m_mcToSelfMap.clear();
    m_caloHitToMCMap.clear();
    m_mcToCaloHitListMap.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingEventShowerStartFinderTool::FillTrueHitMaps(const MCParticleList *const pMCParticleList, const CaloHitList *const pCaloHitList)
{
    LArMCParticleHelper::GetMCToSelfMap(pMCParticleList, m_mcToSelfMap);
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(pCaloHitList, m_mcToSelfMap, m_caloHitToMCMap, m_mcToCaloHitListMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingEventShowerStartFinderTool::GetTargetShowers(const MCParticleList *const pMCParticleList, TargetHierarchyMap &targetHierarchyMap) const
{
    MCParticleVector targetShowerVector;

    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        // Is it an EM particle
        if (!this->IsEM(pMCParticle))
            continue;

        // Does it not exist within an EM hierarchy?
        const MCParticleList &parentMCParticleList(pMCParticle->GetParentList());

        if (parentMCParticleList.size() == 0)
            continue;

        if (parentMCParticleList.size() != 1)
        {
            std::cout << "Parent MCParticle list is more than one??" << std::endl;
            throw;
        }

        if (this->IsEM(parentMCParticleList.front()))
            continue;

        if (m_mcToCaloHitListMap.find(pMCParticle) == m_mcToCaloHitListMap.end())
            continue;

        targetShowerVector.push_back(pMCParticle);
    }

    for (const MCParticle *const pMCParticle : targetShowerVector)
    {
        MCParticleVector particleHierarchy;
        this->FillHierarchyVector(pMCParticle, particleHierarchy);

        if (particleHierarchy.empty())
            continue;

        unsigned int hitCount(m_mcToCaloHitListMap.at(pMCParticle).size());

        for (const MCParticle *const pHierarchyParticle : particleHierarchy)
        {
            if (m_mcToCaloHitListMap.find(pHierarchyParticle) == m_mcToCaloHitListMap.end())
                continue;

            hitCount += m_mcToCaloHitListMap.at(pHierarchyParticle).size();
        }

        if (hitCount < m_threshold2DHitCount)
            continue;

        targetHierarchyMap[pMCParticle] = particleHierarchy;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingEventShowerStartFinderTool::IsEM(const MCParticle *const pMCParticle) const
{
    const int pdgCode(std::abs(pMCParticle->GetParticleId()));

    return ((pdgCode == 11) || (pdgCode == 22));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingEventShowerStartFinderTool::FillHierarchyVector(const MCParticle *const pParentMCParticle, MCParticleVector &particleHierarchy) const
{
    for (const MCParticle *const pMCParticle : pParentMCParticle->GetDaughterList())
    {
        if (std::find(particleHierarchy.begin(), particleHierarchy.end(), pMCParticle) != particleHierarchy.end())
            continue;

        particleHierarchy.push_back(pMCParticle);

        this->FillHierarchyVector(pMCParticle, particleHierarchy);
    }
}

    /*

void CheatingShowerStartFinderTool::GetEMShowerLead(const MCParticle *const pMCParticle, const MCParticle *&pParentMCParticle)
{
    pParentMCParticle = pMCParticle;

    do
    {
        const MCParticleList &parentMCParticleList(pMCParticle->GetParentList());

        if (parentMCParticleList.size() == 0)
            return;

        if (parentMCParticleList.size() != 1)
        {
            std::cout << "Parent MCParticle list is more than one??" << std::endl;
            throw;
        }

        pParentMCParticle = parentMCParticleList.front();
    }
    while (this->IsEM(pParentMCParticle));
}

    */

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingEventShowerStartFinderTool::GetViewShowerStart(const MCParticle *const pParentMCParticle, const MCParticleVector &hierarchyParticles, 
    const HitType hitType, CartesianVector &viewShowerStart) const
{
    const CartesianVector &primaryVertex(pParentMCParticle->GetVertex());
    const CartesianVector viewPrimaryVertex(LArGeometryHelper::ProjectPosition(this->GetPandora(), primaryVertex, hitType));
    const CartesianVector &primaryDirection(pParentMCParticle->GetMomentum().GetUnitVector());
    const CartesianVector viewPrimaryDirection(LArGeometryHelper::ProjectDirection(this->GetPandora(), primaryDirection, hitType));

    // Find closest clear shower branch to primary vertex
    bool found(false);
    float closestSepSq(std::numeric_limits<float>::max());

    for (const MCParticle *const pHierarchyParticle : hierarchyParticles)
    {
        if (m_mcToCaloHitListMap.find(pHierarchyParticle) == m_mcToCaloHitListMap.end())
            continue;

        // Collect branch hits
        unsigned int viewCaloHitCount(0);
        const CaloHitList &caloHitList(m_mcToCaloHitListMap.at(pHierarchyParticle));

        for (const CaloHit *const pCaloHit : caloHitList)
        {
            if (pCaloHit->GetHitType() != hitType)
                continue;

            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

            if (pHierarchyParticle->GetParticleId() == 22)
            {
                const float separationSquared(viewPrimaryDirection.GetCrossProduct(hitPosition - viewPrimaryVertex).GetMagnitudeSquared());

                if (separationSquared < (m_hitDistanceThreshold * m_hitDistanceThreshold))
                    continue;

                ++viewCaloHitCount;
            }
        }

        if (viewCaloHitCount < m_branchHitThreshold)
            continue;

        // Is it the closest?
        const CartesianVector &vertex(pHierarchyParticle->GetVertex());
        const CartesianVector viewVertex(LArGeometryHelper::ProjectPosition(this->GetPandora(), vertex, hitType));
        const float primaryVertexSepSq((viewPrimaryVertex - viewVertex).GetMagnitudeSquared());

        if (primaryVertexSepSq < closestSepSq)
        {
            found = true;
            closestSepSq = primaryVertexSepSq;
            viewShowerStart = viewVertex;
        }
    }

    // Adjust so that vertex lies on shower axis...
    const float longitudinalProjection(viewPrimaryDirection.GetDotProduct(viewShowerStart - viewPrimaryVertex));
    viewShowerStart = (viewPrimaryVertex + (viewPrimaryDirection * longitudinalProjection));

    return found;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingEventShowerStartFinderTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "Threshold2DHitCount", m_threshold2DHitCount));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "HitDistanceThreshold", m_hitDistanceThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "BranchHitThreshold", m_branchHitThreshold));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
