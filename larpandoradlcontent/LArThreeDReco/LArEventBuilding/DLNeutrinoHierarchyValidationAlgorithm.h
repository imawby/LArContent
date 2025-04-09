/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLNeutrinoHierarchyValidationAlgorithm.h
 *
 *  @brief  Header file for the DL neutrino hierarchy validation algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_NEUTRINO_HIERARCHY_VALIDATION_ALGORITHM_H
#define LAR_DL_NEUTRINO_HIERARCHY_VALIDATION_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"

#include "larpandoradlcontent/LArCheating/DLCheatHierarchyTool.h"
#include "larpandoradlcontent/LArThreeDReco/LArEventBuilding/LArHierarchyPfo.h"

using namespace lar_content;

namespace lar_dl_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DLNeutrinoHierarchyValidationAlgorithm class
 */
class DLNeutrinoHierarchyValidationAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DLNeutrinoHierarchyValidationAlgorithm();

    /**
     *  @brief  Default destructor
     */
    ~DLNeutrinoHierarchyValidationAlgorithm();

private:
    typedef std::map<const pandora::ParticleFlowObject *, const pandora::MCParticle *> PfoToMCParticleMap;
    typedef std::map<const pandora::ParticleFlowObject *, std::pair<const pandora::ParticleFlowObject *, int>> ChildToParentPfoMap;

    pandora::StatusCode Run();

    /**
     *  @brief  Return the neutrino pfo
     *
     *  @param  the pointer to the neutrino pfo to fill
     *
     *  @return whether the neutrino pfo can be found
     */
    bool GetNeutrinoPfo(const pandora::ParticleFlowObject *&pNeutrinoPfo) const;

    /**
     *  @brief  Fill the track and shower-like HierarchyPfoVector
     *
     *  @param  pNeutrinoPfo a pointer to the neutrino pfo
     *  @param  trackPfos the track-like HierarchyPfoVector
     *  @param  showerPfos the shower-like HierarchyPfoVector
     */
    void FillTrackShowerVectors(const pandora::ParticleFlowObject *const pNeutrinoPfo, HierarchyPfoVector &trackPfos, HierarchyPfoVector &showerPfos) const;

    float GetNuVertexAccuracy(const pandora::ParticleFlowObject *const pNeutrinoPfo) const;

    void Validate(const pandora::ParticleFlowObject *const pNeutrinoPfo, const HierarchyPfoVector &trackPfos,
                  const HierarchyPfoVector &showerPfos) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_validationFileName;
    std::string m_validationTreeName;   
    std::string m_mcParticleListName;                   ///< the name of the MCParticle list
    std::string m_neutrinoPfoListName;                  ///< the name of the neutrino pfo list
    pandora::StringVector m_pfoListNames;               ///< the name of the pfo lists
    unsigned int m_minClusterSize;                      ///< the minimum threshold of 3D hits of a considered pfo
    int m_slidingFitWindow;                             ///< the sliding fit window to use in pfo sliding linear fits
    DLCheatHierarchyTool *m_cheatHierarchyTool;         ///< The tool used to obtain the true hierarchy
};

} // namespace lar_dl_content

#endif // #ifndef LAR_DL_NEUTRINO_HIERARCHY_VALIDATION_ALGORITHM_H
