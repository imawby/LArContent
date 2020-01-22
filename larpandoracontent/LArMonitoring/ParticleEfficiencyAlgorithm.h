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

    void FillTreeWithSingleParticleEvent(const pandora::MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoCompletenessMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoPurityMap);

    void FillTreeWithTwoParticleEvent(const pandora::MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoCompletenessMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoPurityMap);

    void AddMatchesEntryToTree(const pandora::MCParticle *const pMCParticle, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoCompletenessMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoPurityMap);

    void FillTreeWithUnmatchedSingleParticleEvent(const pandora::MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap);

    void FillTreeWithUnmatchedTwoParticleEvent(const pandora::MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap);

    void AddNoPfoEntryToTree(const pandora::MCParticle *const pMCParticle, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap);

    void AddMCParticleDataToTree(const pandora::MCParticle *const pMCParticle, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap);

    void VisualizeReconstructableMCParticles(const pandora::MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap);

    void VisualizeReconstructedPfos(const pandora::PfoVector &orderedPfoVector, const LArMCParticleHelper::PfoContributionMap &pfoToRecoHitsMap);

    void PrintMCParticleMatchesInfoToScreen(const pandora::MCParticleVector &orderedTargetMCParticleVector, const pandora::PfoVector &orderedPfoVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoCompletenessMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoPurityMap);

    void GetLArSoftAngles(const pandora::CartesianVector &vector, float &theta0XZ, float &theta0YZ);

    void GetDeltaLArSoftAngles(const pandora::MCParticle *const particle1, const pandora::MCParticle *const particle2, float &deltaTheta0XZ, float &deltaTheta0YZ);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string  m_caloHitListName; // Name of input calo hit list
    std::string  m_pfoListName; // Name of input pfo list

    int m_particlesInEvent;

    LArMCParticleHelper::PrimaryParameters m_parameters;
    bool m_foldToPrimaries; ///< whether to fold all hits to primary pfos and MC particles

    float m_minCompleteness;
    float m_minPurity;

    bool m_writeToTree;
    std::string m_treeName;
    std::string m_fileName;
    bool m_printToScreen;
    bool m_visualiseMCParticles;
    bool m_visualisePfos;

    int m_eventNumber;

  };

} //namespace lar_content


#endif // LAR_PFO_VALIDATION_ALGORITHM_H
