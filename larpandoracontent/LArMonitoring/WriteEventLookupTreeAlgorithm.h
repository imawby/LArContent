/**
 *  @file   larpandoracontent/LArMonitoring/WriteEventLookupTreeAlgorithm.h
 *
 *  @brief Header file for the write event lookup tree algorithm
 *
 * $Log: $
 */

#ifndef LAR_WRITE_EVENT_LOOKUP_TREE_ALGORITHM_H
#define LAR_WRITE_EVENT_LOOKUP_TREE_ALGORITHM_H

#include "Pandora/Algorithm.h"


namespace lar_content {

  /**
   * @brief WriteEventLookupTreeAlgorithm class
   */
  class WriteEventLookupTreeAlgorithm : public pandora::Algorithm 
  {

  public:

      WriteEventLookupTreeAlgorithm();

      ~WriteEventLookupTreeAlgorithm();

  private:

      pandora::StatusCode Run();
      void FillTreeWithParticleEvent(const pandora::MCParticleVector &orderedTargetMCParticleVector);
      void AddMCParticleDataToTree(const pandora::MCParticle *const pMCParticle);
      void GetLArSoftAngles(const pandora::CartesianVector &vector, float &theta0XZ, float &theta0YZ);
      pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


      bool m_writeToTree;
      std::string m_treeName;
      std::string m_fileName;
      int m_eventNumber;
  };


} //namespace lar_content


#endif // LAR_WRITE_EVENT_LOOKUP_TREE_ALGORITHM_H
