/**
 *  @file   larpandoracontent/LArCheating/CheatingShowerStartFinderAlgorithm.h
 *
 *  @brief  Header file for the cheating shower start finder class.
 *
 *  $Log: $
 */
#ifndef LAR_CHEATING_SHOWER_START_FINDER_ALGORITHM_H
#define LAR_CHEATING_SHOWER_START_FINDER_ALGORITHM_H 1

#include "larpandoracontent/LArControlFlow/MasterAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  CheatingShowerStartFinderAlgorithm class
 */
class CheatingShowerStartFinderAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingShowerStartFinderAlgorithm();

    pandora::StatusCode Run();

private:
    typedef std::map<int, pandora::MCParticleList> HierarchyMap;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void FillHierarchyVector(const pandora::MCParticle *const pParentMCParticle, pandora::MCParticleVector &particleHierarchy) const;

    void PrepareHierarchyVector(const pandora::MCParticle *const pMCPrimary, 
        const LArMCParticleHelper::MCContributionMap &mcToCaloHitListMap, const pandora::MCParticleVector &particleHierarchy,
        pandora::MCParticleVector &filteredHierarchyU, pandora::MCParticleVector &filteredHierarchyV, pandora::MCParticleVector &filteredHierarchyW) const;

    std::string m_mcParticleListName;
    std::string m_caloHitListName;
    float m_hitDistanceThreshold;
    unsigned int m_branchHitThreshold;
    float m_branchAngleThreshold;
    float m_branchDistanceThreshold;
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_SHOWER_START_FINDER_ALGORITHM_H
