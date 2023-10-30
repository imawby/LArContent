/**
 *  @file   larpandoracontent/LArShowerRefinement/CheatingShowerStartFinderTool.cc
 *
 *  @brief  Implementation of the shower start finder tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArCheating/CheatingShowerStartFinderTool.h"


using namespace pandora;

namespace lar_content
{

CheatingShowerStartFinderTool::CheatingShowerStartFinderTool() : 
    m_mcToCaloHitListMap(),
    m_threshold2DHitCount(100),
    m_hitDistanceThreshold(1.f),
    m_branchHitThreshold(5),
    m_visualize(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
// calo hit list is the shower spine hit list
StatusCode CheatingShowerStartFinderTool::Run(const ParticleFlowObject *const pShowerPfo, const MCParticleList *const pMCParticleList, 
    const CaloHitList &caloHitList, const HitType hitType, CartesianVector &viewShowerStart)
{
    const MCParticle *pTargetMCParticle(nullptr);

    if (this->GetTargetShower(pShowerPfo, hitType, pTargetMCParticle) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;

    /////////////////////////////////////////
    if (m_visualize)
    {
        CaloHitList showerCaloHits;
        LArPfoHelper::GetCaloHits(pShowerPfo, hitType, showerCaloHits);
        PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &showerCaloHits, "showerCaloHits", RED);
        PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &caloHitList, "showerSpine", VIOLET);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }
    /////////////////////////////////////////

    // Do truth information for the shower spine
    this->ResetMaps();
    this->FillTrueHitMaps(pShowerPfo, hitType, pMCParticleList, caloHitList);

    MCParticleVector particleHierarchy;
    this->FillHierarchyVector(pTargetMCParticle, particleHierarchy);

    if (this->GetViewShowerStart(caloHitList, pTargetMCParticle, particleHierarchy, hitType, viewShowerStart))
        return STATUS_CODE_SUCCESS;

    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingShowerStartFinderTool::GetTargetShower(const ParticleFlowObject *const pShowerPfo, const HitType &hitType, 
    const MCParticle *&pTargetMCParticle) const
{
    CaloHitList showerCaloHits;
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, showerCaloHits);

    std::map<const MCParticle*, int> contributionMap;
    int highestHitNumber(-std::numeric_limits<int>::max());
    double highestEnergy(-std::numeric_limits<double>::max());
    const MCParticle *pMainMCParticle(nullptr);

    for (const CaloHit *const pCaloHit : showerCaloHits)
    {
        try
        {
            const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));

            if (this->IsEM(pMCParticle))
            {
                const MCParticle *pLeadMCParticle(pMCParticle);
                this->GetEMShowerLead(pMCParticle, pLeadMCParticle);

                pMCParticle = pLeadMCParticle;
            }

            if (contributionMap.find(pMCParticle) == contributionMap.end())
                contributionMap[pMCParticle] = 1;
            else
                contributionMap[pMCParticle] += 1;


            if ((contributionMap[pMCParticle] > highestHitNumber) || 
               ((contributionMap[pMCParticle] == highestHitNumber) && (pMCParticle->GetEnergy() > highestEnergy)))
            {
                highestHitNumber = contributionMap[pMCParticle];
                highestEnergy = pMCParticle->GetEnergy();
                pMainMCParticle = pMCParticle;
            }
        }
        catch (const StatusCodeException &)
        {
        }
    }

    if (!pMainMCParticle)
        return STATUS_CODE_NOT_FOUND;

    if ((pMainMCParticle->GetParticleId() != 22) && (pMainMCParticle->GetParticleId() != 11))
        return STATUS_CODE_NOT_FOUND;

    pTargetMCParticle = pMainMCParticle;

    return STATUS_CODE_SUCCESS;
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

void CheatingShowerStartFinderTool::ResetMaps()
{
    m_mcToCaloHitListMap.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CheatingShowerStartFinderTool::FillTrueHitMaps(const ParticleFlowObject *const pShowerPfo, const HitType hitType, 
    const MCParticleList *const pMCParticleList, const CaloHitList &showerSpineHitList)
{
    // Combine HitLists
    CaloHitList combinedHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, combinedHitList);

    for (const CaloHit *const pCaloHit : showerSpineHitList)
    {
        if (std::find(combinedHitList.begin(), combinedHitList.end(), pCaloHit) == combinedHitList.end())
            combinedHitList.push_back(pCaloHit);
    }

    // Fill m_mcToCaloHitListMap
    for (const CaloHit *const pCaloHit : combinedHitList)
    {
        try
        {
            const MCParticleWeightMap &weightMap(pCaloHit->GetMCParticleWeightMap());

            for (auto &entry : weightMap)
                m_mcToCaloHitListMap[entry.first].push_back(pCaloHit);
        }
        catch (const StatusCodeException &)
        {
        }
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

bool CheatingShowerStartFinderTool::GetViewShowerStart(const CaloHitList &showerSpine, const MCParticle *const pParentMCParticle, 
    const MCParticleVector &hierarchyParticles, const HitType hitType, CartesianVector &viewShowerStart) const
{
    try
    {
        // Fit the spine...
        CartesianPointVector showerSpinePositions;

        for (const CaloHit *const pJam : showerSpine)
            showerSpinePositions.push_back(pJam->GetPositionVector());

        const TwoDSlidingFitResult spineFitResult(&showerSpinePositions, 20, LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));

        const CartesianVector &primaryVertex(pParentMCParticle->GetVertex());
        const CartesianVector viewPrimaryVertex(LArGeometryHelper::ProjectPosition(this->GetPandora(), primaryVertex, hitType));
        const CartesianVector &primaryDirection(pParentMCParticle->GetMomentum().GetUnitVector());
        const CartesianVector viewPrimaryDirection(LArGeometryHelper::ProjectDirection(this->GetPandora(), primaryDirection, hitType));

        // Find closest clear shower branch to primary vertex
        bool found(false);
        float closestSepSq(std::numeric_limits<float>::max());
        CaloHitList candidateHits;

        for (const MCParticle *const pHierarchyParticle : hierarchyParticles)
        {
            if (m_mcToCaloHitListMap.find(pHierarchyParticle) == m_mcToCaloHitListMap.end())
                continue;

            // Collect branch hits
            const CartesianVector &vertex(pHierarchyParticle->GetVertex());
            const CartesianVector viewVertex(LArGeometryHelper::ProjectPosition(this->GetPandora(), vertex, hitType));
            CaloHitVector offAxisHits;

            for (const CaloHit *const pCaloHit : m_mcToCaloHitListMap.at(pHierarchyParticle))
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
               const float separationSquared(viewPrimaryDirection.GetCrossProduct(hitPosition - viewPrimaryVertex).GetMagnitudeSquared());

               if (separationSquared < (0.5f * 0.5f))
                    continue;

                offAxisHits.push_back(pCaloHit);
            }

            // Need branch to be significant
            if (offAxisHits.size() < m_branchHitThreshold)
                continue;

            // Sort by distance to vertex to find a closest connected branch.. 
            std::sort(offAxisHits.begin(), offAxisHits.end(),
            LArConnectionPathwayHelper::SortByDistanceToPoint(viewVertex));

            // Find closest branch (this is inefficient isobel)
            bool foundBranch(false);
            CartesianVector branchVertex(0.f, 0.f, 0.f);
            CaloHitList branchHits;

            for (const CaloHit *const pSeedCaloHit : offAxisHits)
            {
                CaloHitList branch;
                branch.push_back(pSeedCaloHit);

                bool hitsAdded(true);

                while (hitsAdded)
                {
                    hitsAdded = false;

                    for (const CaloHit *const pCaloHit : offAxisHits)
                    {
                        if (std::find(branch.begin(), branch.end(), pCaloHit) != branch.end())
                            continue;

                        const float closestDistance(LArHitWidthHelper::GetClosestDistance(pCaloHit, branch));

                        if (closestDistance < m_hitDistanceThreshold)
                        {
                            branch.push_back(pCaloHit);
                            hitsAdded = true;
                        }
                    }
                }

                if (branch.size() >= m_branchHitThreshold)
                {
                    branchHits.insert(branchHits.begin(), branch.begin(), branch.end());
                    foundBranch = true;
                    branchVertex = branch.front()->GetPositionVector();

                    break;
                }
            }

            if (!foundBranch)
                continue;

            // Is it the closest?
            const float primaryVertexSepSq((viewPrimaryVertex - branchVertex).GetMagnitudeSquared());

            if (primaryVertexSepSq < closestSepSq)
            {
                found = true;
                closestSepSq = primaryVertexSepSq;
                viewShowerStart = branchVertex;
                candidateHits.clear();
                candidateHits.insert(candidateHits.begin(), branchHits.begin(), branchHits.end());
            }
        }

        if (!found)
            return false;

        // Adjust so that vertex lies on spine axis...
        float rL(0.f), rT(0.f);
        spineFitResult.GetLocalPosition(viewShowerStart, rL, rT);

        ///////////////////////
        if (m_visualize)
        {
            std::cout << "rL: " << rL << std::endl;
            std::cout << "rT: " << rT << std::endl;
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &viewShowerStart, "viewShowerStart", BLUE, 2);
            PandoraMonitoringApi::ViewEvent(this->GetPandora());
        }

        if (spineFitResult.GetGlobalFitPosition(rL, viewShowerStart) != STATUS_CODE_SUCCESS)
            return false;

        int minLayer(spineFitResult.GetMinLayer());
        int maxLayer(spineFitResult.GetMaxLayer());

        if ((rL > spineFitResult.GetL(maxLayer)) || (rL < spineFitResult.GetL(minLayer)))
            return false;

        if (rT > 3.f)
            return false;

        return true;
    }
    catch (const StatusCodeException &)
    {
        return false;
    }
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

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "Visualize", m_visualize));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
