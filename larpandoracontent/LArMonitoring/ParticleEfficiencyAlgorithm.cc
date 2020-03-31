/**
 *  @file   larpandoracontent/LArMonitoring/ParticleEfficiencyAlgorithm
 *
 *  @brief  Implementation of the performance assessment algorithm
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArMonitoring/ParticleEfficiencyAlgorithm.h"



using namespace pandora;

namespace lar_content {

  ParticleEfficiencyAlgorithm::ParticleEfficiencyAlgorithm() :
    m_caloHitListName(), 
    m_pfoListNames(),
    m_vertexListName(),
    m_particlesInEvent(),
    m_foldToPrimaries(false), 
    m_minCompleteness(0.7),
    m_minPurity(0.7),
    m_writeToTree(false), 
    m_treeName(),
    m_fileName(), 
    m_printToScreen(false),
    m_visualiseMCParticles(false),
    m_visualisePfos(false),
    m_isThreeViewMode(true),
    m_eventNumber(0)
  {
  }

//------------------------------------------------------------------------------------------------------------------------------------------

  ParticleEfficiencyAlgorithm::~ParticleEfficiencyAlgorithm() 
  {
    if(m_writeToTree) {
      try {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName, m_fileName, "UPDATE"));
	std::cout << "TREE SAVED" << std::endl;
      }
      catch (const StatusCodeException &) {
        std::cout << "UNABLE TO WRITE TREE" << std::endl; 
      }
    }
  }

//------------------------------------------------------------------------------------------------------------------------------------------


  StatusCode ParticleEfficiencyAlgorithm::Run() {

    ++m_eventNumber;

    std::cout << "Event Number: " << m_eventNumber << std::endl;

    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));
   
    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));


    // set reconstruction parameters to be appropriate for a single view
    if(!m_isThreeViewMode)
    {
        m_parameters.m_minPrimaryGoodHits = 5;
        m_parameters.m_minHitsForGoodView = 5;
        m_parameters.m_minPrimaryGoodViews = 1;
    }

    // MC PARTICLES

    // Get MC Particle to reconstructable hits map
    LArMCParticleHelper::MCContributionMap mcToRecoHitsMap;
    LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, m_parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, mcToRecoHitsMap, m_foldToPrimaries);

    // For user output purposes
    MCParticleVector orderedTargetMCParticleVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({mcToRecoHitsMap}, orderedTargetMCParticleVector);


    // VISUALISE ON SCREEN
    // Visualize the reconstructable MC particles and their hits
    if(m_visualiseMCParticles) {
      std::cout << "MC PARTICLES AND THEIR RECONSTRUCTABLE HITS (ON DISPLAY) " << std::endl;
      VisualizeReconstructableMCParticles(orderedTargetMCParticleVector, mcToRecoHitsMap);
    }

    // PFOS

    // WRITE TO TREE

    // There are two lists, track particles and shower particles 
    // Get all the pfos from each list (depending on which ones exist in Pandora...)
    PfoList allPfos;
    unsigned int listDoesNotExistCount(0);
    for(const std::string &listName : m_pfoListNames) {

      const PfoList *pPfoList = nullptr;
      StatusCode listStatus = PandoraContentApi::GetList(*this, listName, pPfoList);

      if(listStatus == STATUS_CODE_NOT_INITIALIZED) {
        listDoesNotExistCount++;
        continue;
      } 

      // Ensure otther StatusCodes are still handled
      PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, listStatus);
    
      for(const Pfo *const pPfo : *pPfoList)
        allPfos.push_back(pPfo);

    }

    // If no pfos are created, or the lists have not been created yet, still need to fill MC tree with this info
    if((listDoesNotExistCount == m_pfoListNames.size()) || (allPfos.empty())){
      if(m_writeToTree) {
	if(m_particlesInEvent == 1) {
	  FillTreeWithUnmatchedSingleParticleEvent(orderedTargetMCParticleVector, mcToRecoHitsMap);
	}
	if(m_particlesInEvent == 2) {
	  FillTreeWithUnmatchedTwoParticleEvent(orderedTargetMCParticleVector, mcToRecoHitsMap);
	}
      }
      return STATUS_CODE_SUCCESS;
    }
 

    // Get pfo to reconstructable hits map
    LArMCParticleHelper::PfoContributionMap pfoToRecoHitsMap;

    if (m_foldToPrimaries) {
      // Get list of 'primary' pfos to be matched with the target MC particles
      PfoList finalStatePfos;
      for (const ParticleFlowObject *const pPfo : allPfos) {
	      if (LArPfoHelper::IsFinalState(pPfo))
	          finalStatePfos.push_back(pPfo);
      }
      
      LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(finalStatePfos, mcToRecoHitsMap, pfoToRecoHitsMap, true);
    } else {
        LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(allPfos, mcToRecoHitsMap, pfoToRecoHitsMap, false);
    }

    // For user output purposes
    PfoVector orderedPfoVector;
    LArMonitoringHelper::GetOrderedPfoVector(pfoToRecoHitsMap, orderedPfoVector);

    // VISUALISE ON SCREEN
    // Visualize the reconstructed pfos and their hits
    if(m_visualisePfos) {
      std::cout << "RECONSTRUCTED PFOS AND THEIR RECONSTRUCTABLE HITS (ON DISPLAY) " << std::endl;
      VisualizeReconstructedPfos(orderedPfoVector, pfoToRecoHitsMap);
    }


    // MC PARTICLE AND PFO

    // Find hits that they share
    LArMCParticleHelper::PfoToMCParticleHitSharingMap pfoToMCParticleHitSharingMap;
    LArMCParticleHelper::MCParticleToPfoHitSharingMap mcParticleToPfoHitSharingMap;
    LArMCParticleHelper::GetPfoMCParticleHitSharingMaps(pfoToRecoHitsMap, {mcToRecoHitsMap}, pfoToMCParticleHitSharingMap, mcParticleToPfoHitSharingMap);

    // Calculate purity and completeness for MC->Pfo matches
    LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap mcParticleToPfoCompletenessMap;
    LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap mcParticleToPfoPurityMap; 
    LArMCParticleHelper::GetMCToPfoCompletenessPurityMaps(mcToRecoHitsMap, pfoToRecoHitsMap, mcParticleToPfoHitSharingMap, mcParticleToPfoCompletenessMap, mcParticleToPfoPurityMap);

    // WRITE TO TREE
    if(m_writeToTree) {
      if(m_particlesInEvent == 1) {
	FillTreeWithSingleParticleEvent(orderedTargetMCParticleVector, mcToRecoHitsMap, mcParticleToPfoHitSharingMap, mcParticleToPfoCompletenessMap, mcParticleToPfoPurityMap);
      }
      if(m_particlesInEvent == 2) {
	FillTreeWithTwoParticleEvent(orderedTargetMCParticleVector, mcToRecoHitsMap, mcParticleToPfoHitSharingMap, mcParticleToPfoCompletenessMap, mcParticleToPfoPurityMap);
      }
    }

    // PRINT MATCHES INFO TO SCREEN

    if(m_printToScreen) {
      PrintMCParticleMatchesInfoToScreen(orderedTargetMCParticleVector, orderedPfoVector, mcToRecoHitsMap, mcParticleToPfoHitSharingMap, mcParticleToPfoCompletenessMap, mcParticleToPfoPurityMap);
    }
    

    return STATUS_CODE_SUCCESS;

  }


//------------------------------------------------------------------------------------------------------------------------------------------

  void ParticleEfficiencyAlgorithm::FillTreeWithSingleParticleEvent(const MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoCompletenessMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoPurityMap) {

    MCParticleVector targetPrimaries;

    for(const MCParticle *const pMCParticle : orderedTargetMCParticleVector) {
      if(LArMCParticleHelper::IsPrimary(pMCParticle))
	targetPrimaries.push_back(pMCParticle);
    }

    if(targetPrimaries.size() != 1) {
      if(targetPrimaries.size() > 1) {
        std::cout << "WARNING: EVENT HAS MORE THAT 1 PRIMARY PARTICLE" << std::endl;
      }
      
      std::cout << "WARNING: EVENT HAS LESS THAN 1 PRIMARY PARTICLE, RETURNED WIHTOUT WRITING TO TREE" << std::endl;
      return;
    }


    for(const MCParticle *const pMCParticle : targetPrimaries) {
      AddMatchesEntryToTree(pMCParticle, mcToRecoHitsMap, mcParticleToPfoHitSharingMap, mcParticleToPfoCompletenessMap, mcParticleToPfoPurityMap);
    }

  }

//------------------------------------------------------------------------------------------------------------------------------------------ 


  void ParticleEfficiencyAlgorithm::FillTreeWithTwoParticleEvent(const MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoCompletenessMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoPurityMap) {

    MCParticleVector targetPrimaries;

    for(const MCParticle *const pMCParticle : orderedTargetMCParticleVector) {
      if(LArMCParticleHelper::IsPrimary(pMCParticle))
	targetPrimaries.push_back(pMCParticle);
    }
 
    if(targetPrimaries.size() != 2) {
      if(targetPrimaries.size() > 2) {
	std::cout << "WARNING: EVENT HAS MORE THAT 2 PRIMARY PARTICLES" << std::endl;
      }

      std::cout << "WARNING: EVENT HAS LESS THAN 2 PRIMARY PARTICLES, RETURNED WIHTOUT WRITING TO TREE" << std::endl;

      return;
    }

    float openingAngle = targetPrimaries[0]->GetMomentum().GetOpeningAngle(targetPrimaries[1]->GetMomentum());
    
    for(const MCParticle *const pMCParticle : targetPrimaries) {
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "OpeningAngle", openingAngle));
      AddMatchesEntryToTree(pMCParticle, mcToRecoHitsMap, mcParticleToPfoHitSharingMap, mcParticleToPfoCompletenessMap, mcParticleToPfoPurityMap);
    }
    
  }

//------------------------------------------------------------------------------------------------------------------------------------------

  void ParticleEfficiencyAlgorithm::AddMatchesEntryToTree(const MCParticle *const pMCParticle, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoCompletenessMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoPurityMap) {

    AddMCParticleDataToTree(pMCParticle, mcToRecoHitsMap);

    std::vector<int> sharedHitsVector;
    std::vector<double> completenessVector;
    std::vector<double> purityVector;

    // If no matches made
    if(mcParticleToPfoHitSharingMap.at(pMCParticle).empty()) {
      sharedHitsVector.push_back(0);
      completenessVector.push_back(0);
      purityVector.push_back(0);
    }

    for(auto pfoSharedHitPair : mcParticleToPfoHitSharingMap.at(pMCParticle)) {

      sharedHitsVector.push_back(pfoSharedHitPair.second.size());

      for(auto pfoCompletenessPair : mcParticleToPfoCompletenessMap.at(pMCParticle)) {
        if(pfoSharedHitPair.first == pfoCompletenessPair.first)
          completenessVector.push_back(pfoCompletenessPair.second);
      }

      for(auto pfoPurityPair : mcParticleToPfoPurityMap.at(pMCParticle)) {
        if(pfoSharedHitPair.first == pfoPurityPair.first)
          purityVector.push_back(pfoPurityPair.second);
      }
    }

    int isReconstructed(0);
    for(unsigned int i(0); i < completenessVector.size(); ++i) {
      if((completenessVector[i] > m_minCompleteness) && (purityVector[i] > m_minPurity)) {
	isReconstructed = 1;
	break;
      }
    }

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "IsReconstructed", isReconstructed));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "SharedHitsVector", &sharedHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "CompletenessVector", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "PurityVector", &purityVector));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
  
  }

//------------------------------------------------------------------------------------------------------------------------------------------

  void ParticleEfficiencyAlgorithm::FillTreeWithUnmatchedSingleParticleEvent(const MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap) {

    MCParticleVector targetPrimaries;

    for(const MCParticle *const pMCParticle : orderedTargetMCParticleVector) {
      if(LArMCParticleHelper::IsPrimary(pMCParticle))
	targetPrimaries.push_back(pMCParticle);
    }

    if(targetPrimaries.size() != 1) {
      if(targetPrimaries.size() > 1) {
        std::cout << "WARNING: EVENT HAS MORE THAT 1 PRIMARY PARTICLE" << std::endl;
      }
      
      std::cout << "WARNING: EVENT HAS LESS THAN 1 PRIMARY PARTICLE, RETURNED WIHTOUT WRITING TO TREE" << std::endl;
      return;
    }

    for(const MCParticle *const pMCParticle : targetPrimaries) {
        AddNoPfoEntryToTree(pMCParticle, mcToRecoHitsMap);
    }

  }


//------------------------------------------------------------------------------------------------------------------------------------------ 

  void ParticleEfficiencyAlgorithm::FillTreeWithUnmatchedTwoParticleEvent(const MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap) {

    MCParticleVector targetPrimaries;

    for(const MCParticle *const pMCParticle : orderedTargetMCParticleVector) {
      if(LArMCParticleHelper::IsPrimary(pMCParticle))
	targetPrimaries.push_back(pMCParticle);
    }

    if(targetPrimaries.size() != 2) {
      if(targetPrimaries.size() > 2) {
        std::cout << "WARNING: EVENT HAS MORE THAT 2 PRIMARY PARTICLES" << std::endl;
      }
      std::cout << "WARNING: EVENT HAS LESS THAN 2 PRIMARY PARTICLES, RETURNED WIHTOUT WRITING TO TREE" << std::endl;
      return;
    }

    float openingAngle = targetPrimaries[0]->GetMomentum().GetOpeningAngle(targetPrimaries[1]->GetMomentum());

    for(const MCParticle *const pMCParticle : targetPrimaries) {
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "OpeningAngle", openingAngle));
      AddNoPfoEntryToTree(pMCParticle, mcToRecoHitsMap);
    }

  }


//------------------------------------------------------------------------------------------------------------------------------------------ 

  void ParticleEfficiencyAlgorithm::AddNoPfoEntryToTree(const MCParticle *const pMCParticle, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap) {

    AddMCParticleDataToTree(pMCParticle, mcToRecoHitsMap);

    std::vector<int> sharedHitsVector({0});
    std::vector<double> completenessVector({0});
    std::vector<double> purityVector({0});

    int isReconstructed(0);

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "SharedHitsVector", &sharedHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "CompletenessVector", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "PurityVector", &purityVector));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "IsReconstructed", isReconstructed));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
    
  }

//------------------------------------------------------------------------------------------------------------------------------------------

  void ParticleEfficiencyAlgorithm::AddMCParticleDataToTree(const MCParticle *const pMCParticle, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap) {

    float theta0YZ;
    float theta0XZ;
    GetLArSoftAngles(pMCParticle->GetMomentum(), theta0XZ, theta0YZ);


    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "EventNumber", m_eventNumber));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "MCParticleID", pMCParticle->GetParticleId()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Energy", pMCParticle->GetEnergy()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Momentum", pMCParticle->GetMomentum().GetMagnitude()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "X", pMCParticle->GetVertex().GetX()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Y", pMCParticle->GetVertex().GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Z", pMCParticle->GetVertex().GetZ()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Hierarchy", LArMCParticleHelper::GetHierarchyTier(pMCParticle)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Theta0YZ", theta0YZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Theta0XZ", theta0XZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "TotHits", static_cast<int>(mcToRecoHitsMap.at(pMCParticle).size())));


    if(!m_vertexListName.empty())
        GetReconstructedOffsetFromEventVertex(pMCParticle);

  }



//------------------------------------------------------------------------------------------------------------------------------------------

  void ParticleEfficiencyAlgorithm::VisualizeReconstructableMCParticles(const MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap) {

    // Visualize the target MC particles
    PandoraMonitoringApi::Create(this->GetPandora());
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);
    srand(time(NULL));
    for(const MCParticle *const pMCParticle : orderedTargetMCParticleVector) {
      CaloHitList uHits;
      CaloHitList vHits;
      CaloHitList wHits;
      for(const CaloHit *const caloHit : mcToRecoHitsMap.at(pMCParticle)) { 
	if(caloHit->GetHitType() == TPC_VIEW_U) {
	  uHits.push_back(caloHit);
	} else if(caloHit->GetHitType() == TPC_VIEW_V) {
	  vHits.push_back(caloHit);
	} else {
	  wHits.push_back(caloHit);
	}
      }
      Color colour = static_cast<Color>((rand() % (Color::LIGHTYELLOW - 1)) + 1);
      std::string name = "PDG: " + std::to_string(pMCParticle->GetParticleId()) + " Hierarchy Tier: " + std::to_string(LArMCParticleHelper::GetHierarchyTier(pMCParticle));
      PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &uHits, name + " (" + std::to_string(uHits.size()) + " U HITS)", colour);
      PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &vHits, name + " (" + std::to_string(vHits.size()) + " V HITS)", colour);
      PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &wHits, name + " (" + std::to_string(wHits.size()) + " W HITS)", colour);

      PandoraMonitoringApi::Pause(this->GetPandora());
    }

    PandoraMonitoringApi::ViewEvent(this->GetPandora());

  }

//------------------------------------------------------------------------------------------------------------------------------------------


  void ParticleEfficiencyAlgorithm::VisualizeReconstructedPfos(const PfoVector &orderedPfoVector, const LArMCParticleHelper::PfoContributionMap &pfoToRecoHitsMap) {

    typedef std::map<const ParticleFlowObject*, int> PfoToIdMap; 

    PfoToIdMap pfoToIdMap;
    for(unsigned int id(0); id < orderedPfoVector.size(); ++id) {
      pfoToIdMap[orderedPfoVector[id]] = (id + 1);
    }
    
    // Visualize Pfos
    PandoraMonitoringApi::Create(this->GetPandora());
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_DEFAULT, -1.f, 1.f, 1.f);
    srand(time(NULL));

    for(const ParticleFlowObject *const pPfo: orderedPfoVector) {
      CaloHitList uHits;
      CaloHitList vHits;
      CaloHitList wHits;
      for(const CaloHit *const caloHit : pfoToRecoHitsMap.at(pPfo)) { 
	if(caloHit->GetHitType() == TPC_VIEW_U) {
	  uHits.push_back(caloHit);
	} else if(caloHit->GetHitType() == TPC_VIEW_V) {
	  vHits.push_back(caloHit);
	} else {
	  wHits.push_back(caloHit);
	}
      }
      Color colour = static_cast<Color>((rand() % (Color::LIGHTYELLOW - 1)) + 1);
      std::string name = "Id: " + std::to_string(pfoToIdMap.at(pPfo)) + " Hierarchy Tier: " + std::to_string(LArPfoHelper::GetHierarchyTier(pPfo));
      PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &uHits, name + " (" + std::to_string(uHits.size()) + " U HITS)", colour);
      PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &vHits, name + " (" + std::to_string(vHits.size()) + " V HITS)", colour);
      PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &wHits, name + " (" + std::to_string(wHits.size()) + " W HITS)", colour);

      PandoraMonitoringApi::Pause(this->GetPandora());
    }

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
  }



//------------------------------------------------------------------------------------------------------------------------------------------


  void ParticleEfficiencyAlgorithm::PrintMCParticleMatchesInfoToScreen(const MCParticleVector &orderedTargetMCParticleVector, const PfoVector &orderedPfoVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoCompletenessMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoPurityMap) {
  
    // Label pfos
    typedef std::map<const ParticleFlowObject*, int> PfoToIdMap;

    PfoToIdMap pfoToIdMap;
    for(unsigned int id(0); id < orderedPfoVector.size(); ++id) {
      pfoToIdMap[orderedPfoVector[id]] = (id + 1);
    }

    unsigned int reconstructedMCParticles(0);

    for(const MCParticle *const pMCParticle : orderedTargetMCParticleVector) {

      // Get completeness and purity for all particles matched to a MC particle
      LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap::const_iterator completenessIter = mcParticleToPfoCompletenessMap.find(pMCParticle);
      LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap::const_iterator purityIter = mcParticleToPfoPurityMap.find(pMCParticle);

      bool isReconstructed(false);

      int uHits(0);
      int vHits(0);
      int wHits(0);
      for(const CaloHit *const caloHit : mcToRecoHitsMap.at(pMCParticle)) { 
	if(caloHit->GetHitType() == TPC_VIEW_U) {
	  ++uHits;
	} else if(caloHit->GetHitType() == TPC_VIEW_V) {
	  ++vHits;
	} else {
	  ++wHits;
	}
      }

      float theta0XZ;
      float theta0YZ;
      GetLArSoftAngles(pMCParticle->GetMomentum(), theta0XZ, theta0YZ);

      std::string mcString = "MC Particle: (PDG: " + std::to_string(pMCParticle->GetParticleId()) + " Hierarchy Tier: " + std::to_string(LArMCParticleHelper::GetHierarchyTier(pMCParticle)) + " uHits: " + std::to_string(uHits) + " vHits: " + std::to_string(vHits) + " wHits: " + std::to_string(wHits) + ")" + "\n" + "Theta0XZ: " + std::to_string(theta0XZ*180/M_PI) + " Theta0YZ: " + std::to_string(theta0YZ*180/M_PI);
      std::cout << mcString << std::endl;

      std::cout << completenessIter->second.size() << " match(es) made: ";

      if(!completenessIter->second.size()) {
	std::cout << std::endl;
	std::cout << "\033[31m" << "NOT RECONSTRUCTED" << "\033[0m" << std::endl;
	continue;
      }
      
      std::cout << "(Pfo Id, Shared Hits, Completeness, Purity)" << std::endl;

      for(LArMCParticleHelper::PfoToSharedHitsVector::const_iterator iter(mcParticleToPfoHitSharingMap.at(pMCParticle).begin()); iter != mcParticleToPfoHitSharingMap.at(pMCParticle).end(); ++iter) {

	std::cout << "(" << pfoToIdMap.at(iter->first) << ", " << iter->second.size();

	// This step finds the purity and completeness for the same pfo, vectors are looped over incase lists are not in the same order
        LArMCParticleHelper::PfoCompletenessPurityPair matchedPfoCompletenessPair;
        LArMCParticleHelper::PfoCompletenessPurityPair matchedPfoPurityPair;
        for(LArMCParticleHelper::PfoToCompletenessPurityVector::const_iterator pfoCompletenessPairIter(completenessIter->second.begin()); pfoCompletenessPairIter != completenessIter->second.end(); ++pfoCompletenessPairIter) {
	  if(pfoToIdMap.at(pfoCompletenessPairIter->first) == pfoToIdMap.at(iter->first))   
	    matchedPfoCompletenessPair = *pfoCompletenessPairIter;
	}
        for(LArMCParticleHelper::PfoToCompletenessPurityVector::const_iterator pfoPurityPairIter(purityIter->second.begin()); pfoPurityPairIter != purityIter->second.end(); ++pfoPurityPairIter) {
	  if(pfoToIdMap.at(pfoPurityPairIter->first) == pfoToIdMap.at(iter->first)) 
	    matchedPfoPurityPair = *pfoPurityPairIter;
	}
	std::cout << ", " << matchedPfoCompletenessPair.second;
	std::cout << ", " << matchedPfoPurityPair.second << ")" << std::endl;

	if((matchedPfoCompletenessPair.second > 0.7) && (matchedPfoPurityPair.second > 0.7)) {
	  isReconstructed = true;
	  reconstructedMCParticles++;
	}
      
      }

      isReconstructed ? std::cout <<  "\033[32m" << "RECONSTRUCTED" << "\033[0m" << std::endl : std::cout << "\033[31m" << "NOT RECONSTRUCTED" << "\033[0m" << std::endl;
	
    }  

      std::cout << "Reconstruction Efficiency: " << static_cast<double>(reconstructedMCParticles*100)/static_cast<double>(orderedTargetMCParticleVector.size()) << "%" << std::endl;

  }

//------------------------------------------------------------------------------------------------------------------------------------------


  void ParticleEfficiencyAlgorithm::GetLArSoftAngles(const CartesianVector &vector, float &theta0XZ, float &theta0YZ) {

    theta0YZ = asin(vector.GetY()/vector.GetMagnitude());
    theta0XZ = atan2(vector.GetX(), vector.GetZ());

  }

//------------------------------------------------------------------------------------------------------------------------------------------



void ParticleEfficiencyAlgorithm::GetReconstructedOffsetFromEventVertex(const MCParticle *const pMCParticle) 
{

    const VertexList *pVertexList = nullptr;
    StatusCode vertexStatus = PandoraContentApi::GetList(*this, m_vertexListName, pVertexList);

    float vertexOffset(0);

    // If vertex is not reconstructed, still want to be able to fill the tree
    // Fill will something non-sensical so can retrieve this information
    if(vertexStatus != STATUS_CODE_SUCCESS)
    {
        if(vertexStatus == STATUS_CODE_NOT_INITIALIZED)
        {
            vertexOffset = -1.f;
            PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "RecoOffsetFromEventVertex", vertexOffset));
	    return;
        }
        else
        {
            // Ensure that other SatusCodes are still handled
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, vertexStatus);
        }
    }

    // Check that there is not more than one reconstructed event vertex (this shouldn't happen)
    if(pVertexList->size() != 1)
    {
        std::cout << "WARNING - VERTEX LIST SIZE DOES NOT EQUAL 1" << std::endl;
        throw;
    }

    // Reconstructed event vertex
    const Vertex *const pRecoNeutrinoVertex = (*pVertexList->begin());
    const CartesianVector *const pRecoEventVertex(&pRecoNeutrinoVertex->GetPosition()); 

    // MC particle vertex
    const CartesianVector *const pVertexPosition(&pMCParticle->GetVertex());

    vertexOffset = std::sqrt(pVertexPosition->GetDistanceSquared(*pRecoEventVertex));

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "RecoOffsetFromEventVertex", vertexOffset));

}

//------------------------------------------------------------------------------------------------------------------------------------------

  void ParticleEfficiencyAlgorithm::GetDeltaLArSoftAngles(const MCParticle *const particle1, const MCParticle *const particle2, float &deltaTheta0XZ, float &deltaTheta0YZ) {

    float particle1Theta0XZ;
    float particle1Theta0YZ;

    float particle2Theta0XZ;
    float particle2Theta0YZ;

    GetLArSoftAngles(particle1->GetMomentum(), particle1Theta0XZ, particle1Theta0YZ);
    GetLArSoftAngles(particle2->GetMomentum(), particle2Theta0XZ, particle2Theta0YZ);

    deltaTheta0XZ = fabs(particle1Theta0XZ - particle2Theta0XZ);
    deltaTheta0YZ = fabs(particle1Theta0YZ - particle2Theta0YZ);

  }


//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

  StatusCode ParticleEfficiencyAlgorithm::ReadSettings(const TiXmlHandle xmlHandle) {

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadVectorOfValues(xmlHandle,
        "PfoListNames", m_pfoListNames));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "VertexListName", m_vertexListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ParticlesInEvent", m_particlesInEvent));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FoldToPrimaries", m_foldToPrimaries));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCompleteness", m_minCompleteness));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinPurity", m_minPurity));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "TreeName", m_treeName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "FileName", m_fileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintToScreen", m_printToScreen));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VisualiseMCParticles", m_visualiseMCParticles));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "VisualisePfos", m_visualisePfos));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "IsThreeViewMode", m_isThreeViewMode));



  return STATUS_CODE_SUCCESS;

  }

} //namespace lar_content

