#include "Pandora/AlgorithmHeaders.h"


#include "larpandoracontent/LArMonitoring/PullDataAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArInteractionTypeHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

using namespace pandora;

namespace lar_content {

  PullDataAlgorithm::PullDataAlgorithm() :
    m_eventNumber(0),
    m_writeToTree(false),
    m_treeName(),
    m_treeFileName(),
    m_writeToFile(false),
    m_fileName(),
    m_PDG()
  {  
  }

  PullDataAlgorithm::~PullDataAlgorithm() {
    if(m_writeToTree) {
      try {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName, m_treeFileName, "UPDATE"));
      }
      catch (const StatusCodeException &) {
        std::cout << "UNABLE TO WRITE TREE" << std::endl;
      }
    }
  }

  StatusCode PullDataAlgorithm::Run() {

    m_eventNumber++;
    std::cout << "Event Number: "  << m_eventNumber << std::endl;

    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    const CaloHitList *pCaloHitList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));



    // Interaction type for the event 
    MCParticleVector primaryParticlesVector;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaryParticlesVector);
    MCParticleList primaryParticlesList;
    std::copy(primaryParticlesVector.begin(), primaryParticlesVector.end(), std::back_inserter(primaryParticlesList));
    LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(primaryParticlesList));

    std::cout << "Interaction: " << LArInteractionTypeHelper::ToString(interactionType) << std::endl;

    /*
    // Get MC Particle to reconstructable hits map
    LArMCParticleHelper::MCContributionMap mcToRecoHitsMap;
    LArMCParticleHelper::SelectUnfoldedReconstructableMCParticles(pMCParticleList, pCaloHitList, m_parameters, mcToRecoHitsMap);
    
    
    if(m_writeToTree) {
        WriteToParticleEventTree(pMCParticleList, mcToRecoHitsMap);
      //WriteToMuonProtonEventTree(pMCParticleList);
    }  

    if(m_writeToFile) {
      //WriteToMuonProtonEventFile(pMCParticleList);
    }
    */
   return STATUS_CODE_SUCCESS;

  }

//------------------------------------------------------------------------------------------------------------------------------------------

  void PullDataAlgorithm::GetLArSoftAngles(const CartesianVector &vector, float &theta0XZ, float &theta0YZ) {

    theta0YZ = asin(vector.GetY()/vector.GetMagnitude());
    theta0XZ = atan2(vector.GetX(), vector.GetZ());

  }

//------------------------------------------------------------------------------------------------------------------------------------------

  bool PullDataAlgorithm::IsParticleReconstructable(const MCParticle *const pMCParticle, LArMCParticleHelper::MCContributionMap mcToRecoHitsMap) {

    LArMCParticleHelper::MCContributionMap::iterator iter = mcToRecoHitsMap.find(pMCParticle);
    if(iter == mcToRecoHitsMap.end()) {
      return false; 
    }

      return true;
  }

//------------------------------------------------------------------------------------------------------------------------------------------

  void PullDataAlgorithm::WriteToParticleEventTree(const MCParticleList *pMCParticleList, LArMCParticleHelper::MCContributionMap mcToRecoHitsMap) {

    // Nuance Code for event
    const MCParticle *pMCNeutrino(nullptr);
    /*
    for(const MCParticle *const pMCParticle : *pMCParticleList) {
        
      if(LArMCParticleHelper::IsDownstreamOfBeamNeutrino(pMCParticle)) {
	pMCNeutrino = LArMCParticleHelper::GetParentMCParticle(pMCParticle);
	break;
      }
    }
    */
    unsigned int nuanceCode(LArMCParticleHelper::GetNuanceCode(pMCNeutrino));

    // Interaction type for the event 
    MCParticleVector primaryParticlesVector;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaryParticlesVector);
    MCParticleList primaryParticlesList;
    std::copy(primaryParticlesVector.begin(), primaryParticlesVector.end(), std::back_inserter(primaryParticlesList));
    LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(primaryParticlesList));


    for(const MCParticle *const pMCParticle : *pMCParticleList) {

      // For Muon distribution determination
      if(abs(pMCParticle->GetParticleId()) != m_PDG)
        continue;

      if(!IsParticleReconstructable(pMCParticle, mcToRecoHitsMap)) 
	continue;

      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "InteractionType", static_cast<int>(interactionType)));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "NuanceCode", static_cast<int>(nuanceCode)));

      WriteMCParticleToTree(pMCParticle);

      PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
    }

    return;

  }

//------------------------------------------------------------------------------------------------------------------------------------------

  
  void PullDataAlgorithm::WriteToMuonProtonEventTree(const MCParticleList *pMCParticleList) {

    // Count how many primary muons and protons are in interaction
    unsigned int numMuons(0);
    unsigned int numProtons(0);
    MCParticleVector targetPrimaries;

    for(const MCParticle *pMCParticle : *pMCParticleList) {

      // Check whether primary and part of beam neutrino interaction
      if(!LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCParticle))
	continue;

      if(abs(pMCParticle->GetParticleId()) == 13) {
	numMuons++;
	targetPrimaries.push_back(pMCParticle);
      }

      if(abs(pMCParticle->GetParticleId()) == 2212) {
	numProtons++;
	targetPrimaries.push_back(pMCParticle);
      }

    }

    // If not an event with a single primary proton and muon then return
    if(!((numMuons == 1) && (numProtons == 1)))
      return;

    // Nuance Code for event
    const MCParticle *pMCNeutrino(nullptr);
    /*
    for(const MCParticle *const pMCParticle : *pMCParticleList) {
      if(LArMCParticleHelper::IsDownstreamOfBeamNeutrino(pMCParticle)) {
	pMCNeutrino = LArMCParticleHelper::GetParentMCParticle(pMCParticle);
	break;
      }
    }
    */
    unsigned int nuanceCode(LArMCParticleHelper::GetNuanceCode(pMCNeutrino));

    // Interaction type for the event 
    MCParticleVector primaryParticlesVector;
    LArMCParticleHelper::GetPrimaryMCParticleList(pMCParticleList, primaryParticlesVector);
    MCParticleList primaryParticlesList;
    std::copy(primaryParticlesVector.begin(), primaryParticlesVector.end(), std::back_inserter(primaryParticlesList));
    LArInteractionTypeHelper::InteractionType interactionType(LArInteractionTypeHelper::GetInteractionType(primaryParticlesList));


    float angleBetween;
    angleBetween = targetPrimaries[0]->GetMomentum().GetOpeningAngle(targetPrimaries[1]->GetMomentum());

    for(const MCParticle *const pMCParticle : targetPrimaries) {

      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "InteractionType", static_cast<int>(interactionType)));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "NuanceCode", static_cast<int>(nuanceCode)));
      PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "OpeningAngle", angleBetween));

      WriteMCParticleToTree(pMCParticle);

      /*
      float theta0UZ;
      float theta0VZ;
      this->GetTheta0UZ(pMCParticle->GetMomentum(), theta0UZ);
      this->GetTheta0VZ(pMCParticle->GetMomentum(), theta0VZ);

      std::cout << "MC Particle: " << pMCParticle->GetParticleId() << std::endl;
      std::cout << "Theta0UZ: " << theta0UZ*180/M_PI << std::endl;
      std::cout << "Theta0VZ: " << theta0VZ*180/M_PI << std::endl;
      */
      PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
    }

    return;

  }
  

//------------------------------------------------------------------------------------------------------------------------------------------

  void PullDataAlgorithm::WriteMCParticleToTree(const MCParticle *const pMCParticle) {
    
    float theta0XZ;
    float theta0YZ;

    GetLArSoftAngles(pMCParticle->GetMomentum(), theta0XZ, theta0YZ);

    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "EventNumber", m_eventNumber));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "ParticleID", pMCParticle->GetParticleId()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Momentum", pMCParticle->GetMomentum().GetMagnitude()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Hierarchy", LArMCParticleHelper::GetHierarchyTier(pMCParticle)));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Energy", pMCParticle->GetEnergy()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "X", pMCParticle->GetVertex().GetX()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Y", pMCParticle->GetVertex().GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Z", pMCParticle->GetVertex().GetZ()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "PX", pMCParticle->GetMomentum().GetX()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "PY", pMCParticle->GetMomentum().GetY()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "PZ", pMCParticle->GetMomentum().GetZ()));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Theta0XZ", theta0XZ));
    PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), m_treeName, "Theta0YZ", theta0YZ));

    return;
  } 


//------------------------------------------------------------------------------------------------------------------------------------------

  void PullDataAlgorithm::WriteToMuonProtonEventFile(const MCParticleList *pMCParticleList) {

    // Count how many primary muons and protons are in interaction
    unsigned int numMuons(0);
    unsigned int numProtons(0);
    MCParticleVector targetPrimaries;

    for(const MCParticle *pMCParticle : *pMCParticleList) {

      // Check whether primary and part of beam neutrino interaction
      if(!LArMCParticleHelper::IsBeamNeutrinoFinalState(pMCParticle))
	continue;

      if(abs(pMCParticle->GetParticleId()) == 13) {
	numMuons++;
	targetPrimaries.push_back(pMCParticle);
      }

      if(abs(pMCParticle->GetParticleId()) == 2212) {
	numProtons++;
	targetPrimaries.push_back(pMCParticle);
      }
    }

    // If not an event with a single primary proton and muon then return
    if(!((numMuons == 1) && (numProtons == 1)))
      return;

    std::ofstream outfile(m_fileName, std::ios::app);
    outfile << m_eventNumber << " 2" << std::endl;

    for(const MCParticle *pMCParticle : targetPrimaries) {
      CartesianVector momentum(pMCParticle->GetMomentum());
      CartesianVector vertex(pMCParticle->GetVertex());
      outfile << "1 " << pMCParticle->GetParticleId() << " 0 0 0 0 " << momentum.GetX() << " " << momentum.GetY() << " " << momentum.GetZ() << " " << pMCParticle->GetEnergy();
      pMCParticle->GetParticleId() == 13 ? outfile << " 0.105 " : outfile << " 0.938 ";
      outfile << vertex.GetX()/10 << " " << vertex.GetY()/10 << " " << vertex.GetZ()/10 << " 0" << std::endl;
    }

    outfile.close();

  }

//------------------------------------------------------------------------------------------------------------------------------------------


 StatusCode PullDataAlgorithm::ReadSettings(const TiXmlHandle xmlHandle) {

   PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeName", m_treeName));
   PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteToTree", m_writeToTree));
   PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "TreeFileName", m_treeFileName));
   PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "WriteToFile", m_writeToFile));
   PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FileName", m_fileName));
   PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "PDG", m_PDG));
   return STATUS_CODE_SUCCESS;

 }

} //namespace lar_content


/**

    for(const MCParticle *const pMCParticle : *pMCParticleList) {

      if(!((pMCParticle->GetParticleId() == 13) || (abs(pMCParticle->GetParticleId()) == 2212)))
	continue;

      std::cout << "ParticleID: " << pMCParticle->GetParticleId() << std::endl;
      std::cout << "PX: " << pMCParticle->GetMomentum().GetX() << std::endl;
      std::cout << "PY: " << pMCParticle->GetMomentum().GetY() << std::endl;
      std::cout << "PZ: " << pMCParticle->GetMomentum().GetZ() << std::endl;
      std::cout << "Energy: " << pMCParticle->GetEnergy() << std::endl;
      std::cout << "Mass: " << std::endl;
      std::cout << "X: " << pMCParticle->GetVertex().GetX() << std::endl;
      std::cout << "Y: " << pMCParticle->GetVertex().GetY() << std::endl;
      std::cout << "Z: " << pMCParticle->GetVertex().GetZ() << std::endl;

    }
 */
