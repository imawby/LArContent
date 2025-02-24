/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/TrackShowerScoreValidationAlgorithm.h
 *
 *  @brief  Header file for the track shower score validation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_MLP_TRACK_SHOWER_SCORE_VALIDATION_ALGORITHM_H
#define LAR_MLP_TRACK_SHOWER_SCORE_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoradlcontent/LArCheating/MLPCheatHierarchyTool.h"

using namespace lar_content;

namespace lar_dl_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  TrackShowerScoreValidationAlgorithm class
 */
class TrackShowerScoreValidationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    TrackShowerScoreValidationAlgorithm();

    /**
     *  @brief  Default destructor
     */
    ~TrackShowerScoreValidationAlgorithm();

private:
    typedef std::map<const pandora::ParticleFlowObject*, const pandora::MCParticle*> PfoToMCParticleMap;
    typedef std::map<const pandora::ParticleFlowObject*, std::pair<const pandora::ParticleFlowObject*, int>> ChildToParentPfoMap;

    pandora::StatusCode Run();

    float GetNSpacepoints(const pandora::ParticleFlowObject *const pPfo) const;    

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    int m_eventID;
    pandora::StringVector m_pfoListNames;    
    std::string m_validationFileName;
    std::string m_validationTreeName;    
    MLPCheatHierarchyTool *m_cheatHierarchyTool;
};

} // namespace lar_dl_content

#endif // #ifndef LAR_MLP_TRACK_SHOWER_SCORE_VALIDATION_ALGORITHM_H
