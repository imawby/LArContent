/**
 *  @file   larpandoracontent/LArShowerRefinement/EventMLShowerInitialRegionRefinementAlgorithm.h
 *
 *  @brief  Header file for the ML shower initial region refinement algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_EVENT_ML_SHOWER_INITIAL_REGION_REFINEMENT_ALGORITHM_H
#define LAR_EVENT_ML_SHOWER_INITIAL_REGION_REFINEMENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArCheating/CheatingEventShowerStartFinderTool.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"

namespace lar_content
{

/**
 *  @brief  EventMLShowerInitialRegionRefinementAlgorithm class
 */
class EventMLShowerInitialRegionRefinementAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    EventMLShowerInitialRegionRefinementAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    void GetShowerStarts(pandora::CartesianPointVector &showerStartsU, pandora::CartesianPointVector &showerStartsV, 
        pandora::CartesianPointVector &showerStartsW) const;

    /**
     *  @brief  Obtain a sorted vector of the reconstructed shower pfos
     *
     *  @param  showerPfoVector the sorted shower vector to be filled
     */
    void FillShowerPfoVector(pandora::PfoVector &showerPfoVector) const;

    /**
     *  @brief  Find and evaluate shower connection pathway, removing if necessary
     *
     *  @param  pShowerPfo the input shower pfo
     */
    void RefineShower(const pandora::ParticleFlowObject *const pShowerPfo) const;

    /**
     *  @brief  Obtain the reconstructed neutrino vertex
     *
     *  @param  nuVertex3D the output neutrino vertex
     *
     *  @return whether a reconstructed neutrino vertex could be found
     */
    pandora::StatusCode GetNeutrinoVertex(pandora::CartesianVector &nuVertex3D) const;

    CheatingEventShowerStartFinderTool *m_pCheatingEventShowerStartFinderTool;           ///< The shower start finder tool

    std::string m_mcParticleListName;
    std::string m_caloHitListName;
    std::string m_showerPfoListName;                           ///< The shower pfo list name
    std::string m_neutrinoVertexListName;                      ///< The neutrino vertex list name
    unsigned int m_minShowerHits3D;                            ///< The min. number of hits of a significant shower
};

} // namespace lar_content

#endif // #ifndef LAR_ML_SHOWER_INITIAL_REGION_REFINEMENT
