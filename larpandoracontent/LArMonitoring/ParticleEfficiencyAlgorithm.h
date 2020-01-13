/**
 *  @file   larpandoracontent/LArMonitoring/ParticleEfficiencyAlgorithm.h
 *
 *  @brief Header file for the particle efficiency algorithm
 *
 * $Log: $
 */

#ifndef LAR_PARTICLE_EFFICIENCY_ALGORITHM_H
#define LAR_PARTICLE_EFFICIENCY_ALGORITHM_H

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content {

  /**
   * @brief ParticleEfficiencyAlgorithm class
   */
  class ParticleEfficiencyAlgorithm : public pandora::Algorithm {

  public:

    /**
     *  @brief Default constructor
     */
    ParticleEfficiencyAlgorithm();

    /**
     *  @brief Default destructor
     */
    ~ParticleEfficiencyAlgorithm();

  private:

    pandora::StatusCode Run();

    void FillTreeWithEvent(const pandora::MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoCompletenessMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoPurityMap);

    void AddMatchesEntryToTree(const pandora::MCParticle *const pMCParticle, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoCompletenessMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoPurityMap);

    void AddNoPfoEntryToTree(const pandora::MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap);

    void AddMCParticleDataToTree(const pandora::MCParticle *const pMCParticle, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap);

    void FillTreeWithTwoParticleEvent(const pandora::MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoCompletenessMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoPurityMap);

    void VisualizeReconstructableMCParticles(const pandora::MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap);

    void VisualizeReconstructedPfos(const pandora::PfoVector &orderedPfoVector, const LArMCParticleHelper::PfoContributionMap &pfoToRecoHitsMap);

    void GetLArSoftAngles(const pandora::CartesianVector &vector, float &theta0XZ, float &theta0YZ);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string  m_caloHitListName; // Name of input calo hit list
    std::string  m_pfoListName; // Name of input pfo list

    LArMCParticleHelper::PrimaryParameters m_parameters;
    bool m_foldToPrimaries; ///< whether to fold all hits to primary pfos and MC particles

    bool m_writeToTree;
    std::string m_treeName;
    std::string m_fileName;
    bool m_printToScreen;

    int m_eventNumber;

  };

} //namespace lar_content


#endif // LAR_PFO_VALIDATION_ALGORITHM_H
