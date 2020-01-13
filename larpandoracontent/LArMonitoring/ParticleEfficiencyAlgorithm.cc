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
    m_pfoListName(),
    m_foldToPrimaries(false), 
    m_writeToTree(false), 
    m_treeName(),
    m_fileName(), 
    m_printToScreen(false),
    m_eventNumber(0)
  {
  }

//------------------------------------------------------------------------------------------------------------------------------------------

  ParticleEfficiencyAlgorithm::~ParticleEfficiencyAlgorithm() 
  {
    if(m_writeToTree) {
      try {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName, m_fileName, "UPDATE"));
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

    // MC PARTICLES

    // Get MC Particle to reconstructable hits map
    LArMCParticleHelper::MCContributionMap mcToRecoHitsMap;
    m_foldToPrimaries ? LArMCParticleHelper::SelectReconstructableMCParticles(pMCParticleList, pCaloHitList, m_parameters, LArMCParticleHelper::IsBeamNeutrinoFinalState, mcToRecoHitsMap) : LArMCParticleHelper::SelectUnfoldedReconstructableMCParticles(pMCParticleList, pCaloHitList, m_parameters, mcToRecoHitsMap);

    // For user output purposes
    MCParticleVector orderedTargetMCParticleVector;
    LArMonitoringHelper::GetOrderedMCParticleVector({mcToRecoHitsMap}, orderedTargetMCParticleVector);

    // PRINT TO SCREEN
    // Visualize the reconstructable MC particles and their hits
    if(m_printToScreen) {
      std::cout << "MC PARTICLES AND THEIR RECONSTRUCTABLE HITS (ON DISPLAY) " << std::endl;
      VisualizeReconstructableMCParticles(orderedTargetMCParticleVector, mcToRecoHitsMap);
    }

    // PFOS

    // WRITE TO TREE
    // If no pfos are created, still need to fill MC tree with this info
    const PfoList *pPfoList = nullptr;
    StatusCode pfoListStatus(PandoraContentApi::GetList(*this, m_pfoListName, pPfoList));

    // Handle case where no pfos are created
    if(pfoListStatus == STATUS_CODE_NOT_INITIALIZED){
      if(m_writeToTree) {
        AddNoPfoEntryToTree(orderedTargetMCParticleVector, mcToRecoHitsMap);
      }
      return STATUS_CODE_NOT_INITIALIZED;
    } 
 
    // ensure that other StatusCodes are still handled
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, pfoListStatus);


    // Get pfo to reconstructable hits map
    LArMCParticleHelper::PfoContributionMap pfoToRecoHitsMap;

    if (m_foldToPrimaries) {
      // Get list of 'primary' pfos to be matched with the target MC particles
      PfoList finalStatePfos;
      for (const ParticleFlowObject *const pPfo : *pPfoList) {
	if (LArPfoHelper::IsFinalState(pPfo))
	  finalStatePfos.push_back(pPfo);
      }
      LArMCParticleHelper::GetPfoToReconstructable2DHitsMap(finalStatePfos, mcToRecoHitsMap, pfoToRecoHitsMap);
    } else {
      LArMCParticleHelper::GetUnfoldedPfoToReconstructable2DHitsMap(*pPfoList, mcToRecoHitsMap, pfoToRecoHitsMap);
    }

    // For user output purposes
    PfoVector orderedPfoVector;
    LArMonitoringHelper::GetOrderedPfoVector(pfoToRecoHitsMap, orderedPfoVector);

    // PRINT TO SCREEN
    if(m_printToScreen) {
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
      //FillTreeWithEvent(orderedTargetMCParticleVector, mcToRecoHitsMap, mcParticleToPfoHitSharingMap, mcParticleToPfoCompletenessMap, mcParticleToPfoPurityMap);
      FillTreeWithTwoParticleEvent(orderedTargetMCParticleVector, mcToRecoHitsMap, mcParticleToPfoHitSharingMap, mcParticleToPfoCompletenessMap, mcParticleToPfoPurityMap);
    }

 
    return STATUS_CODE_SUCCESS;

  }


//------------------------------------------------------------------------------------------------------------------------------------------

  void ParticleEfficiencyAlgorithm::FillTreeWithEvent(const MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap, const LArMCParticleHelper::MCParticleToPfoHitSharingMap &mcParticleToPfoHitSharingMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoCompletenessMap, const LArMCParticleHelper::MCParticleToPfoCompletenessPurityMap &mcParticleToPfoPurityMap) {

    for(const MCParticle *const pMCParticle : orderedTargetMCParticleVector) {
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
 
    if(targetPrimaries.size() != 2)
      return;

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

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "SharedHitsVector", &sharedHitsVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "CompletenessVector", &completenessVector));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "PurityVector", &purityVector));

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
  
  }

//------------------------------------------------------------------------------------------------------------------------------------------

  void ParticleEfficiencyAlgorithm::AddNoPfoEntryToTree(const MCParticleVector &orderedTargetMCParticleVector, const LArMCParticleHelper::MCContributionMap &mcToRecoHitsMap) {
    for(const MCParticle *const pMCParticle : orderedTargetMCParticleVector) {

      AddMCParticleDataToTree(pMCParticle, mcToRecoHitsMap);

      std::vector<int> sharedHitsVector({0});
      std::vector<double> completenessVector({0});
      std::vector<double> purityVector({0});

      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "SharedHitsVector", &sharedHitsVector));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "CompletenessVector", &completenessVector));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "PurityVector", &purityVector));

      PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
    }
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
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Hierarchy", LArMCParticleHelper::GetHierarchyTier(pMCParticle)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Theta0YZ", theta0YZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Theta0XZ", theta0XZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "TotHits", static_cast<int>(mcToRecoHitsMap.at(pMCParticle).size())));

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

  void ParticleEfficiencyAlgorithm::GetLArSoftAngles(const CartesianVector &vector, float &theta0XZ, float &theta0YZ) {

    theta0YZ = asin(vector.GetY()/vector.GetMagnitude());
    theta0XZ = atan2(vector.GetX(), vector.GetZ());

  }

//------------------------------------------------------------------------------------------------------------------------------------------


//------------------------------------------------------------------------------------------------------------------------------------------

  StatusCode ParticleEfficiencyAlgorithm::ReadSettings(const TiXmlHandle xmlHandle) {

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PfoListName", m_pfoListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FoldToPrimaries", m_foldToPrimaries));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "TreeName", m_treeName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "FileName", m_fileName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "PrintToScreen", m_printToScreen));

  return STATUS_CODE_SUCCESS;

  }

} //namespace lar_content

