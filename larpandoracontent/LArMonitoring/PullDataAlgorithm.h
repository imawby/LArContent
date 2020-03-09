/**
 *  @file   larpandoracontent/LArMonitoring/PullDataAlgorithm.h
 *
 *  @brief Header file for the performance assessment algorithm
 *
 * $Log: $
 */

#ifndef LAR_PULL_DATA_ALGORITHM_H
#define LAR_PULL_DATA_ALGORITHM_H

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content {

  class PullDataAlgorithm : public pandora::Algorithm {

  public:

    PullDataAlgorithm();
    ~PullDataAlgorithm();

  private:

    pandora::StatusCode Run();
    void GetLArSoftAngles(const pandora::CartesianVector &vector, float &theta0XZ, float &theta0YZ);
    bool IsParticleReconstructable(const pandora::MCParticle *const pMCParticle, LArMCParticleHelper::MCContributionMap mcToRecoHitsMap);
    void WriteToParticleEventTree(const pandora::MCParticleList *pMCParticleList, LArMCParticleHelper::MCContributionMap mcToRecoHitsMap);
    void WriteToMuonProtonEventTree(const pandora::MCParticleList *pMCParticleList);
    void WriteMCParticleToTree(const pandora::MCParticle *const pMCParticle);
    void WriteToMuonProtonEventFile(const pandora::MCParticleList *pMCParticleList);
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);


    int m_eventNumber;

    bool m_writeToTree;
    std::string m_treeName;
    std::string m_treeFileName;

    bool m_writeToFile;
    std::string m_fileName;

    int m_PDG;

    LArMCParticleHelper::PrimaryParameters m_parameters;

  };



} //namespace lar_content


#endif 