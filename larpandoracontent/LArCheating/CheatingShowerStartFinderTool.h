/**
 *  @file   larpandoracontent/LArShowerRefinement/CheatingShowerStartFinderTool.h
 *
 *  @brief  Header file for the cheating shower start finder tool class.
  *
 *  $Log: $
 */
#ifndef LAR_CHEATING_SHOWER_START_FINDER_TOOL_H
#define LAR_CHEATING_SHOWER_START_FINDER_TOOL_H 1

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

class CheatingShowerStartFinderTool : public pandora::AlgorithmTool
{
public:
    /**
     *  @brief  Default constructor
     */
    CheatingShowerStartFinderTool();

    pandora::StatusCode Run(const pandora::MCParticleList *const pMCParticleList, const pandora::CaloHitList &caloHitList, const pandora::HitType hitType,
        pandora::CartesianVector &viewShowerStart);

private:
    typedef std::map<const pandora::MCParticle*, pandora::MCParticleVector> TargetHierarchyMap;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void ResetMaps();
    void FillTrueHitMaps(const pandora::MCParticleList *const pMCParticleList, const pandora::CaloHitList &caloHitList);
    bool GetTargetShower(const pandora::MCParticle *&pTargetMCParticle) const;
    bool IsEM(const pandora::MCParticle *const pMCParticle) const;
    void GetEMShowerLead(const pandora::MCParticle *const pMCParticle, const pandora::MCParticle *&pParentMCParticle) const;
    void FillHierarchyVector(const pandora::MCParticle *const pParentMCParticle, pandora::MCParticleVector &particleHierarchy) const;
    bool GetViewShowerStart(const pandora::MCParticle *const pParentMCParticle, const pandora::MCParticleVector &hierarchyParticles, 
        const pandora::HitType hitType, pandora::CartesianVector &viewShowerStart) const;

    LArMCParticleHelper::MCRelationMap m_mcToSelfMap;
    LArMCParticleHelper::CaloHitToMCMap m_caloHitToMCMap;
    LArMCParticleHelper::MCContributionMap m_mcToCaloHitListMap; 

    unsigned int m_threshold2DHitCount;
    float m_hitDistanceThreshold;
    unsigned int m_branchHitThreshold;
};

//------------------------------------------------------------------------------------------------------------------------------------------
} // namespace lar_content

#endif // #ifndef LAR_CHEATING_SHOWER_START_FINDER_TOOL_H
