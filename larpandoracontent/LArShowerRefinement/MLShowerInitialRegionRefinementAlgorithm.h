/**
 *  @file   larpandoracontent/LArShowerRefinement/MLShowerInitialRegionRefinementAlgorithm.h
 *
 *  @brief  Header file for the ML shower initial region refinement algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_ML_SHOWER_INITIAL_REGION_REFINEMENT_ALGORITHM_H
#define LAR_ML_SHOWER_INITIAL_REGION_REFINEMENT_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArCheating/CheatingShowerStartFinderTool.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"
#include "larpandoracontent/LArShowerRefinement/PeakDirectionFinderTool.h"
#include "larpandoracontent/LArShowerRefinement/ShowerSpineFinderTool.h"
#include "larpandoracontent/LArShowerRefinement/ShowerStartFinderTool.h"

namespace lar_content
{

/**
 *  @brief  MLShowerInitialRegionRefinementAlgorithm class
 */
class MLShowerInitialRegionRefinementAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MLShowerInitialRegionRefinementAlgorithm();

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
     *  @brief  Build the 2D ProtoShower objects for a given view
     *
     *  @param  pShowerPfo the input shower pfo
     *  @param  nuVertex3D the 3D neutrino vertex
     *  @param  hitType the 2D view
     *  @param  protoShowerVector the output vector of ProtoShower objects
     */
    void BuildViewProtoShowers(const pandora::ParticleFlowObject *const pShowerPfo, const pandora::CartesianVector &nuVertex3D,
        pandora::HitType hitType, ProtoShowerVector &protoShowerVector) const;

    bool GetShowerStart(const pandora::ParticleFlowObject *const pPfo, const pandora::HitType hitType, 
        const pandora::CaloHitList &showerSpineHitList, pandora::CartesianVector &showerStart) const;

    /**
     *  @brief  Obtain the reconstructed neutrino vertex
     *
     *  @param  nuVertex3D the output neutrino vertex
     *
     *  @return whether a reconstructed neutrino vertex could be found
     */
    pandora::StatusCode GetNeutrinoVertex(pandora::CartesianVector &nuVertex3D) const;

    /**
     *  @brief  Obtain the event hit list of a given view
     *
     *  @param  hitType the 2D view
     *  @param  pCaloHitList the output 2D hit list
     *
     *  @return whether a valid 2D hit list could be found
     */
    pandora::StatusCode GetHitListOfType(const pandora::HitType hitType, const pandora::CaloHitList *&pCaloHitList) const;

    PeakDirectionFinderTool *m_pShowerPeakDirectionFinderTool; ///< The shower initial pathway direction finder tool
    ShowerSpineFinderTool *m_pShowerSpineFinderTool;           ///< The shower spine finder tool for the shower
    CheatingShowerStartFinderTool *m_pCheatingShowerStartFinderTool;           ///< The shower start finder tool
    ShowerStartFinderTool *m_pShowerStartFinderTool;           ///< The shower start finder tool

    std::string m_mcParticleListName;
    std::string m_showerPfoListName;                           ///< The shower pfo list name
    std::string m_neutrinoVertexListName;                      ///< The neutrino vertex list name
    std::string m_caloHitListNameU;                            ///< The U calo hit list name
    std::string m_caloHitListNameV;                            ///< The V calo hit list name
    std::string m_caloHitListNameW;                            ///< The W calo hit list name
    unsigned int m_minShowerHits3D;                            ///< The min. number of hits of a significant shower
};

} // namespace lar_content

#endif // #ifndef LAR_ML_SHOWER_INITIAL_REGION_REFINEMENT
