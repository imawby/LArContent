/**
 *  @file   larpandoracontent/LArHelpers/LArMuonLeadingHelper.h
 *
 *  @brief  Header file for the lar muon leading helper class.
 *
 *  $Log: $
 */
#ifndef LAR_MUON_LEADING_HELPER_H
#define LAR_MUON_LEADING_HELPER_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  LArMuonLeadingHelper class
 */
class LArMuonLeadingHelper
{
public:
    /**
     *  @brief   MuonLeadingParameters class
     */
    class ValidationParameters : public LArMCParticleHelper::PrimaryParameters
    {
    public:
        /**
         *  @brief  Constructor
         */
        ValidationParameters();
    };
    
    /**
     *  @brief  Return true if passed a DR tagged MCParticle 
     */    
    static bool IsDeltaRay(const pandora::MCParticle *const pMCParticle);

    static bool IsMichel(const pandora::MCParticle *const pMCParticle);

    static bool IsLeading(const pandora::MCParticle *const pMCParticle);   
    
    static const pandora::MCParticle *GetLeadingParticle(const pandora::MCParticle *const pMCParticle);

    static void GetMCToLeadingMap(const pandora::MCParticleList *const pMCParticleList, LArMCParticleHelper::MCRelationMap &mcToLeadingMap);

    static void RemoveMuonPfosFromList(const pandora::PfoList *const pPfoList, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &unfoldedInterepretedMatchingMap,
        pandora::PfoList &outputList);

    static void SelectNonMuonLeadingPfos(const pandora::PfoList &inputPfoList, pandora::PfoList &outputList);

    static void SelectReconstructableLeadingParticles(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList, const ValidationParameters &parameters,
        LArMCParticleHelper::MCContributionMap &selectedMCParticlesToHitsMap);

    static void GetPfoMatchContamination(const pandora::MCParticle *const pLeadingParticle, const pandora::CaloHitList &matchedPfoHitList,
        pandora::CaloHitList &parentTrackHits, pandora::CaloHitList &otherTrackHits, pandora::CaloHitList &otherShowerHits);

    static void GetMuonPfoContaminationContribution(const pandora::CaloHitList &cosmicRayPfoHitList, const pandora::CaloHitList &leadingMCHitList,
        pandora::CaloHitList &leadingHitsInParentCosmicRay);
    
 private:
    static void SelectLeadingMCParticles(const pandora::MCParticleList *pMCParticleList, pandora::MCParticleVector &selectedParticles);

};

} // namespace lar_content

#endif // #ifndef LAR_MUON_LEADING_HELPER_H
