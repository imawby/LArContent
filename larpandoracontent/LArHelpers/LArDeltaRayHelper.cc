/**
 *  @file   larpandoracontent/LArHelpers/LArDeltaRayHelper.cc
 *
 *  @brief  Implementation of the lar delta ray helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArDeltaRayHelper.h"

#include "Objects/ParticleFlowObject.h"

namespace lar_content
{

using namespace pandora;

LArDeltaRayHelper::DeltaRayParameters::DeltaRayParameters() :
    LArMCParticleHelper::PrimaryParameters(),
    m_maximumContributingTier(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------    

bool LArDeltaRayHelper::IsDeltaRay(const MCParticle *const pMCParticle)
{
    const LArMCParticle *const pLArMCParticle(dynamic_cast<const LArMCParticle*>(pMCParticle));

    // If not ionisation from muon
    if (!pLArMCParticle->GetIsDR())
        return false;
    /*
    const MCParticle *const pParentMuon(LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle));

    // Sanity check
    if (std::fabs(pParentMuon->GetParticleId()) != 13)
        return false;

    const CartesianVector &muonEndpoint(pParentMuon->GetEndpoint());
    const CartesianVector &muonVertex(pParentMuon->GetVertex());
    const CartesianVector &mcParticleVertex(pMCParticle->GetVertex());

    // If parent muon is short
    if ((muonEndpoint - muonVertex).GetMagnitude() < 10.f)
        return false;

    // Reject michel electrons
    if (((muonEndpoint - mcParticleVertex).GetMagnitude() < 2.f) || ((muonVertex - mcParticleVertex).GetMagnitude() < 2.f))
        return false;
    */
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArDeltaRayHelper::IsDecendentOfDeltaRay(const MCParticle *const pMCParticle)
{
    // Navigate upward through MC daughter/parent links - collect this particle and all its parents
    MCParticleVector mcParticleVector;

    const MCParticle *pParentMCParticle = pMCParticle;
    mcParticleVector.push_back(pParentMCParticle);

    while (!pParentMCParticle->GetParentList().empty())
    {
        if (1 != pParentMCParticle->GetParentList().size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pParentMCParticle = *(pParentMCParticle->GetParentList().begin());
        mcParticleVector.push_back(pParentMCParticle);
    }

    // Navigate downward through MC parent/daughter links - search for a delta ray 
    for (MCParticleVector::const_reverse_iterator iter = mcParticleVector.rbegin(), iterEnd = mcParticleVector.rend(); iter != iterEnd; ++iter)
    {
        const MCParticle *const pNextParticle = *iter;

        if (LArDeltaRayHelper::IsDeltaRay(pNextParticle))
            return true;
    }

    return false;
}

    
//------------------------------------------------------------------------------------------------------------------------------------------    

const MCParticle *LArDeltaRayHelper::GetLeadingDeltaRay(const MCParticle *const pMCParticle)
{
    // Navigate upward through MC daughter/parent links - collect this particle and all its parents
    MCParticleVector mcParticleVector;

    const MCParticle *pParentMCParticle = pMCParticle;
    mcParticleVector.push_back(pParentMCParticle);

    while (!pParentMCParticle->GetParentList().empty())
    {
        if (1 != pParentMCParticle->GetParentList().size())
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pParentMCParticle = *(pParentMCParticle->GetParentList().begin());
        mcParticleVector.push_back(pParentMCParticle);
    }

    // Navigate downward through MC parent/daughter links - return the first delta ray 
    for (MCParticleVector::const_reverse_iterator iter = mcParticleVector.rbegin(), iterEnd = mcParticleVector.rend(); iter != iterEnd; ++iter)
    {
        const MCParticle *const pNextParticle = *iter;

        if (LArDeltaRayHelper::IsDeltaRay(pNextParticle))
            return pNextParticle;
    }

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void LArDeltaRayHelper::GetMCToLeadingDeltaRayMap(const MCParticleList *const pMCParticleList, const DeltaRayParameters &parameters,
    LArMCParticleHelper::MCRelationMap &mcToLeadingDeltaRayMap)
{
    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
       if (!LArDeltaRayHelper::IsDecendentOfDeltaRay(pMCParticle))
           continue;

       if (LArMCParticleHelper::GetHierarchyTier(pMCParticle) > parameters.m_maximumContributingTier)
           continue;
       
        mcToLeadingDeltaRayMap[pMCParticle] = LArDeltaRayHelper::GetLeadingDeltaRay(pMCParticle);
    }
}    

//------------------------------------------------------------------------------------------------------------------------------------------

void LArDeltaRayHelper::RemoveMuonPfosFromList(const PfoList *const pPfoList, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &unfoldedInterepretedMatchingMap,
    PfoList &outputList)
{

    std::cout << "JANET" << std::endl;
    
    MCParticleVector reconstructableMCParticles;
    for (auto &entry : unfoldedInterepretedMatchingMap)
        reconstructableMCParticles.push_back(entry.first);

    std::sort(reconstructableMCParticles.begin(), reconstructableMCParticles.end(), LArMCParticleHelper::SortByMomentum);

    PfoList muonPfos;
    for (const MCParticle *const pMCParticle : reconstructableMCParticles)
    {
        if (std::fabs(pMCParticle->GetParticleId()) != 13)
            continue;

        const auto iter(unfoldedInterepretedMatchingMap.find(pMCParticle));

        if (iter == unfoldedInterepretedMatchingMap.end())
            continue;

        std::cout << "BLAH" << std::endl;
        if (!iter->second.empty())
            muonPfos.push_back(iter->second.front().first);
        std::cout << "BLAH 2" << std::endl;
    }

    for (const ParticleFlowObject *const pPfo : *pPfoList)
    {
        if (std::find(muonPfos.begin(), muonPfos.end(), pPfo) == muonPfos.end())
            outputList.push_back(pPfo);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArDeltaRayHelper::SelectNonMuonLeadingPfos(const PfoList &inputPfoList, PfoList &outputList)
{    
    for (const ParticleFlowObject *const pPfo : inputPfoList)
    {
        const PfoList &parentPfoList(pPfo->GetParentPfoList());

        if (parentPfoList.size() > 1)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        if (parentPfoList.empty())
        {
            outputList.push_back(pPfo); 
        }
        else
        {
            if (std::find(inputPfoList.begin(), inputPfoList.end(), parentPfoList.front()) == inputPfoList.end())
                outputList.push_back(pPfo);
        }
    }
}    

//------------------------------------------------------------------------------------------------------------------------------------------    

void LArDeltaRayHelper::SelectReconstructableDeltaRays(const MCParticleList *pMCParticleList, const CaloHitList *pCaloHitList, const DeltaRayParameters &parameters,
    LArMCParticleHelper::MCContributionMap &selectedMCParticlesToHitsMap)
{
    // Obtain map: [mc particle -> target mc particle]
    // ISOBEL - IT DOESN'T MAKE SENSE TO UNFOLD THE HIERARCHY
    LArMCParticleHelper::MCRelationMap mcToTargetMCMap;
    parameters.m_foldBackHierarchy ? LArDeltaRayHelper::GetMCToLeadingDeltaRayMap(pMCParticleList, parameters, mcToTargetMCMap) :
        LArMCParticleHelper::GetMCToSelfMap(pMCParticleList, mcToTargetMCMap);

    // Select reconstructable hits, e.g. remove those downstream of a neutron
    // Unless selectInputHits == false
    // ISOBEL - NOT SURE ABOUT THE PHOTON SHIZZLE
    CaloHitList selectedCaloHitList;
    LArMCParticleHelper::SelectCaloHits(pCaloHitList, mcToTargetMCMap, selectedCaloHitList, parameters.m_selectInputHits, parameters.m_maxPhotonPropagation);

    // Obtain maps: [hit -> target mc particle], [target mc particle -> list of hits]
    LArMCParticleHelper::CaloHitToMCMap trueHitToTargetMCMap;
    LArMCParticleHelper::MCContributionMap targetMCToTrueHitListMap;
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(&selectedCaloHitList, mcToTargetMCMap, trueHitToTargetMCMap, targetMCToTrueHitListMap);

    // Obtain vector: target mc particles
    MCParticleVector targetMCVector;
    LArDeltaRayHelper::SelectTargetMCParticles(pMCParticleList, targetMCVector, parameters);
    
    // Ensure the MCParticles have enough "good" hits to be reconstructed
    LArMCParticleHelper::SelectParticlesByHitCount(targetMCVector, targetMCToTrueHitListMap, mcToTargetMCMap, parameters, selectedMCParticlesToHitsMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArDeltaRayHelper::SelectTargetMCParticles(const MCParticleList *pMCParticleList, MCParticleVector &selectedParticles, const DeltaRayParameters &parameters)
{
    if (parameters.m_foldBackHierarchy)
    {
        for (const MCParticle *const pMCParticle : *pMCParticleList)
        {
            if (LArDeltaRayHelper::IsDeltaRay(pMCParticle))
                selectedParticles.push_back(pMCParticle);
        }
    }
    else
    {
        for (const MCParticle *const pMCParticle : *pMCParticleList)
        {
            if (!LArDeltaRayHelper::IsDecendentOfDeltaRay(pMCParticle))
                continue;
            
            if (LArMCParticleHelper::GetHierarchyTier(pMCParticle) > parameters.m_maximumContributingTier)
                continue;
            
            selectedParticles.push_back(pMCParticle);
        }
    }

    std::sort(selectedParticles.begin(), selectedParticles.end(), LArMCParticleHelper::SortByMomentum);
}

} // namespace lar_content
