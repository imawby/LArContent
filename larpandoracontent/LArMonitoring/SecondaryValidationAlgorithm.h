/**
 *  @file   larpandoracontent/LArMonitoring/SecondaryValidationAlgorithm.h
 *
 *  @brief  Header file for the secondary validation algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_SECONDARY_VALIDATION_ALGORITHM_H
#define LAR_SECONDARY_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  SecondaryValidationAlgorithm class
 */
class SecondaryValidationAlgorithm : public pandora::Algorithm
{
public:
    /**
   *  @brief  Default constructor
   */
    SecondaryValidationAlgorithm();

    virtual ~SecondaryValidationAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void FoldMCNodes(const LArHierarchyHelper::MCHierarchy::Node *const pRootMCNode, LArHierarchyHelper::MCHierarchy::NodeVector &nodeVector);

    void FillTree(const LArHierarchyHelper::MCHierarchy::Node *const pParentMCNode, const pandora::MCParticle *const pMCParent,
                  const LArHierarchyHelper::RecoHierarchy::Node *const pParentRecoNode,
                  const LArHierarchyHelper::MCHierarchy::NodeVector &hierarchyNodes, const LArHierarchyHelper::MCMatchesVector mcMatches,
                  const int tierToExamine);
    
    void FillNullEntry(const pandora::MCParticle *const pMCParent, const pandora::MCParticle *const pMCChild, const int hierarchyTier);

    void FillEntry(const LArHierarchyHelper::MCHierarchy::Node *const pParentMCNode, const LArHierarchyHelper::MCHierarchy::Node *const pChildMCNode, 
                   const LArHierarchyHelper::MCMatches *const pChildMatch, const LArHierarchyHelper::MCMatchesVector &matchesVector, const int hierarchyTier);

    int GetHitsInUpstreamHierarchy(const LArHierarchyHelper::MCHierarchy::Node *const pRootMCNode, 
                                   const LArHierarchyHelper::MCMatchesVector &matchesVector);

    bool FindMCNode(const pandora::MCParticle *const pMCParticle, const LArHierarchyHelper::MCHierarchy::NodeVector &hierarchyMCNodes, 
                    const LArHierarchyHelper::MCHierarchy::Node *&pMCNode);

    bool FindMatch(const pandora::MCParticle *const pMCParticle, const LArHierarchyHelper::MCMatchesVector &matchesVector, 
                   const LArHierarchyHelper::MCMatches *&pRecoMatch);

//void FillTree(const pandora::MCParticle *const pMCParticle, const LArHierarchyHelper::MCHierarchy::NodeVector &mcNodes, const int tierToExamine);
    
    std::string m_caloHitListName;     ///< Name of input calo hit list
    std::string m_pfoListName;         ///< Name of input PFO list
    bool m_writeFile;
    std::string m_fileName;
    std::string m_treeName;
    float m_minPurity;                 ///< Minimum purity to tag a node as being of good quality
    float m_minCompleteness;           ///< Minimum completeness to tag a node as being of good quality
    unsigned int m_minRecoHits;        ///< Minimum number of reconstructed primary good hits
    unsigned int m_minRecoHitsPerView; ///< Minimum number of reconstructed hits for a good view
    unsigned int m_minRecoGoodViews;   ///< Minimum number of reconstructed primary good views
    bool m_removeRecoNeutrons;         ///< Whether to remove reconstructed neutrons and their downstream particles
    bool m_selectRecoHits;             ///< Whether to select reco hits that overlap with the MC particle hits
    int m_maxTierToExamine;
};

} // namespace lar_content

#endif // LAR_SECONDARY_VALIDATION_ALGORITHM_H
