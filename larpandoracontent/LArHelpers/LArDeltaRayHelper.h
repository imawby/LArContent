/**
 *  @file   larpandoracontent/LArHelpers/LArDeltaRayHelper.h
 *
 *  @brief  Header file for the lar delta ray helper class.
 *
 *  $Log: $
 */
#ifndef LAR_DELTA_RAY_HELPER_H
#define LAR_DELTA_RAY_HELPER_H 1

#include "Pandora/PandoraInternal.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  LArDeltaRayHelper class
 */
class LArDeltaRayHelper
{
public:
    /**
     *  @brief   DeltaRayParameters class
     */
    class DeltaRayParameters : public LArMCParticleHelper::PrimaryParameters
    {
    public:
        /**
         *  @brief  Constructor
         */
        DeltaRayParameters();

        unsigned int  m_maximumContributingTier;  
    };
    
    /**
     *  @brief  Return true if passed a DR tagged MCParticle 
     */    
    static bool IsDeltaRay(const pandora::MCParticle *const pMCParticle);

    static bool IsDecendentOfDeltaRay(const pandora::MCParticle *const pMCParticle);    

    static const pandora::MCParticle *GetLeadingDeltaRay(const pandora::MCParticle *const pMCParticle);

    static void GetMCToLeadingDeltaRayMap(const pandora::MCParticleList *const pMCParticleList, const DeltaRayParameters &parameters, LArMCParticleHelper::MCRelationMap &mcToLeadingDeltaRayMap);

    static void RemoveMuonPfosFromList(const pandora::PfoList *const pPfoList, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &unfoldedInterepretedMatchingMap,
        pandora::PfoList &outputList);

    static void SelectNonMuonLeadingPfos(const pandora::PfoList &inputPfoList, pandora::PfoList &outputList);

    static void SelectReconstructableDeltaRays(const pandora::MCParticleList *pMCParticleList, const pandora::CaloHitList *pCaloHitList, const DeltaRayParameters &parameters,
        LArMCParticleHelper::MCContributionMap &selectedMCParticlesToHitsMap);
    
 private:
    static void SelectTargetMCParticles(const pandora::MCParticleList *pMCParticleList, pandora::MCParticleVector &selectedParticles, const DeltaRayParameters &parameters);

};

} // namespace lar_content

#endif // #ifndef LAR_DELTA_RAY_HELPER_H
