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

    void FillHierarchyMap(const pandora::MCParticle *const pParentMCParticle, const int generation, 
        HierarchyMap &hierarchyMap);

    std::string m_mcParticleListName;
    std::string m_caloHitListName;
};

} // namespace lar_content

#endif // #ifndef LAR_CHEATING_SHOWER_START_FINDER_ALGORITHM_H
