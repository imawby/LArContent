/**
 *  @file   larpandoracontent/LArMonitoring/WriteEventLookupTreeAlgorithm
 *
 *  @brief  Implementation of the write event lookup tree algorithm
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"

#include "larpandoracontent/LArMonitoring/WriteEventLookupTreeAlgorithm.h"


using namespace pandora;

namespace lar_content {

WriteEventLookupTreeAlgorithm::WriteEventLookupTreeAlgorithm() :
    m_writeToTree(false), 
    m_treeName("EventLookupTree"),
    m_fileName("EventLookup.root"), 
    m_eventNumber(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

WriteEventLookupTreeAlgorithm::~WriteEventLookupTreeAlgorithm() 
{
    if(m_writeToTree) 
    {
        try 
        {
            PANDORA_MONITORING_API(SaveTree(this->GetPandora(), m_treeName, m_fileName, "UPDATE"));
	    std::cout << "TREE SAVED" << std::endl;
        }
        catch (const StatusCodeException &) 
        {
            std::cout << "UNABLE TO WRITE TREE" << std::endl; 
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------


StatusCode WriteEventLookupTreeAlgorithm::Run() {

    ++m_eventNumber;

    std::cout << "Event Number: " << m_eventNumber << std::endl;

    const MCParticleList *pMCParticleList = nullptr;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pMCParticleList));

    MCParticleVector orderedTargetMCParticleVector(pMCParticleList->begin(), pMCParticleList->end());

    if(m_writeToTree) 
    {
  	FillTreeWithParticleEvent(orderedTargetMCParticleVector);
        
    }

    return STATUS_CODE_SUCCESS;

}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteEventLookupTreeAlgorithm::FillTreeWithParticleEvent(const MCParticleVector &orderedTargetMCParticleVector) {

    MCParticleVector targetPrimaries;

    for(const MCParticle *const pMCParticle : orderedTargetMCParticleVector) 
    {
        if(LArMCParticleHelper::IsPrimary(pMCParticle))
	    targetPrimaries.push_back(pMCParticle);
    }

    if(targetPrimaries.size() < 1) 
    {
        std::cout << "WARNING: EVENT HAS LESS THAN 1 PRIMARY PARTICLE, RETURNED WITHOUT WRITING TO TREE" << std::endl;
        return;
    }


    for(const MCParticle *const pMCParticle : targetPrimaries) 
    {
        AddMCParticleDataToTree(pMCParticle);
    }


}

//------------------------------------------------------------------------------------------------------------------------------------------

void WriteEventLookupTreeAlgorithm::AddMCParticleDataToTree(const MCParticle *const pMCParticle) 
{

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

    PANDORA_MONITORING_API(FillTree(this->GetPandora(), m_treeName));
    
}

//------------------------------------------------------------------------------------------------------------------------------------------


void WriteEventLookupTreeAlgorithm::GetLArSoftAngles(const CartesianVector &vector, float &theta0XZ, float &theta0YZ) 
{

    theta0YZ = asin(vector.GetY()/vector.GetMagnitude());
    theta0XZ = atan2(vector.GetX(), vector.GetZ());

}


//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode WriteEventLookupTreeAlgorithm::ReadSettings(const TiXmlHandle xmlHandle) 
{

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "WriteToTree", m_writeToTree));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TreeName", m_treeName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FileName", m_fileName));

  return STATUS_CODE_SUCCESS;

  }

} //namespace lar_content
