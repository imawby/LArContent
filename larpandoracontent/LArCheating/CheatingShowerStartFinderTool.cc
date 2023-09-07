/**
 *  @file   larpandoracontent/LArShowerRefinement/CheatingShowerStartFinderTool.cc
 *
 *  @brief  Implementation of the shower start finder tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"


#include "larpandoracontent/LArCheating/CheatingShowerStartFinderTool.h"

using namespace pandora;

namespace lar_content
{

CheatingShowerStartFinderTool::CheatingShowerStartFinderTool() : 
    m_threshold2DHitCount(100),
    m_hitDistanceThreshold(1.f),
    m_branchHitThreshold(5)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingShowerStartFinderTool::Run(const MCParticleList *const pMCParticleList, const CaloHitList &caloHitList, const HitType hitType,
    CartesianVector &viewShowerStart)
{
    this->ResetMaps();
    this->FillTrueHitMaps(pMCParticleList, caloHitList);

    const MCParticle *pTargetMCParticle(nullptr);

    if (!this->GetTargetShower(pTargetMCParticle))
        return STATUS_CODE_NOT_FOUND;

    MCParticleVector particleHierarchy;
    this->FillHierarchyVector(pTargetMCParticle, particleHierarchy);

    if (this->GetViewShowerStart(pTargetMCParticle, particleHierarchy, hitType, viewShowerStart))
        return STATUS_CODE_SUCCESS;

    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingShowerStartFinderTool::ResetMaps()
{
    m_mcToSelfMap.clear();
    m_caloHitToMCMap.clear();
    m_mcToCaloHitListMap.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingShowerStartFinderTool::FillTrueHitMaps(const MCParticleList *const pMCParticleList, const CaloHitList &caloHitList)
{
    LArMCParticleHelper::GetMCToSelfMap(pMCParticleList, m_mcToSelfMap);
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(&caloHitList, m_mcToSelfMap, m_caloHitToMCMap, m_mcToCaloHitListMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingShowerStartFinderTool::GetTargetShower(const MCParticle *&pTargetMCParticle) const
{
    bool found(false);

    std::map<const MCParticle*, int> contributionMap;
    int highestHitNumber(-std::numeric_limits<int>::max());
    double highestEnergy(-std::numeric_limits<double>::max());

    for (auto &entry : m_mcToCaloHitListMap)
    {
        const MCParticle *const pMCParticle(entry.first);

        if (!this->IsEM(pMCParticle))
            continue;

        found = true;

        const MCParticle *pEMShowerLead(pMCParticle);
        this->GetEMShowerLead(pMCParticle, pEMShowerLead);

        if (contributionMap.find(pEMShowerLead) == contributionMap.end())
            contributionMap[pEMShowerLead] = entry.second.size();
        else
            contributionMap[pEMShowerLead] += entry.second.size();

        if ((contributionMap[pEMShowerLead] > highestHitNumber) || 
            ((contributionMap[pEMShowerLead] == highestHitNumber) && (pEMShowerLead->GetEnergy() > highestEnergy)))
        {
            highestHitNumber = contributionMap[pEMShowerLead];
            highestEnergy = pEMShowerLead->GetEnergy();
            pTargetMCParticle = pEMShowerLead;
        }
    }

    return found;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CheatingShowerStartFinderTool::IsEM(const MCParticle *const pMCParticle) const
{
    const int pdgCode(std::abs(pMCParticle->GetParticleId()));

    return ((pdgCode == 11) || (pdgCode == 22));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingShowerStartFinderTool::GetEMShowerLead(const MCParticle *const pMCParticle, const MCParticle *&pParentMCParticle) const
{
    if (!this->IsEM(pMCParticle))
        throw;

    const MCParticle *pConsideredMCParticle(pMCParticle);

    while (this->IsEM(pConsideredMCParticle))
    {
        pParentMCParticle = pConsideredMCParticle;

        const MCParticleList &parentMCParticleList(pConsideredMCParticle->GetParentList());

        if (parentMCParticleList.size() == 0)
            return;

        if (parentMCParticleList.size() != 1)
        {
            std::cout << "Parent MCParticle list is more than one??" << std::endl;
            throw;
        }

        pConsideredMCParticle = parentMCParticleList.front();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingShowerStartFinderTool::FillHierarchyVector(const MCParticle *const pParentMCParticle, MCParticleVector &particleHierarchy) const
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

bool CheatingShowerStartFinderTool::GetViewShowerStart(const MCParticle *const pParentMCParticle, const MCParticleVector &hierarchyParticles, 
    const HitType hitType, CartesianVector &viewShowerStart) const
{
    const CartesianVector &primaryVertex(pParentMCParticle->GetVertex());
    const CartesianVector viewPrimaryVertex(LArGeometryHelper::ProjectPosition(this->GetPandora(), primaryVertex, hitType));
    const CartesianVector &primaryDirection(pParentMCParticle->GetMomentum().GetUnitVector());
    const CartesianVector viewPrimaryDirection(LArGeometryHelper::ProjectDirection(this->GetPandora(), primaryDirection, hitType));

    const CartesianVector end(viewPrimaryVertex + (viewPrimaryDirection * 100));
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &viewPrimaryVertex, "viewPrimaryVertex", BLUE, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &viewPrimaryVertex, &end, "viewPrimaryDirection", BLUE, 2, 2);




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
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

            //if (pParentMCParticle->GetParticleId() == 22)
            //{
                const float separationSquared(viewPrimaryDirection.GetCrossProduct(hitPosition - viewPrimaryVertex).GetMagnitudeSquared());

                if (separationSquared < (m_hitDistanceThreshold * m_hitDistanceThreshold))
                    continue;
                //}

            ++viewCaloHitCount;
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

StatusCode CheatingShowerStartFinderTool::ReadSettings(const TiXmlHandle xmlHandle)
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
