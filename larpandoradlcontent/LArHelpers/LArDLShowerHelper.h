/**
 *  @file   larpandoradlcontent/LArHelpers/LArDLShowerHelper.h
 *
 *  @brief  Header file for the lar deep learning helper helper class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_SHOWER_HELPER_H
#define LAR_DL_SHOWER_HELPER_H 1

#include "Pandora/StatusCodes.h"

#include "Objects/CaloHit.h"

#include <set>

namespace lar_dl_content
{

/**
 *  @brief  LArDLShowerHelper class
 */
class LArDLShowerHelper
{
public:
    struct HitFeatures
    {
        HitFeatures();

        float m_xRel;
        float m_zRel;
        float m_rRel;
        float m_cosThetaRel;
        float m_sinThetaRel;
        float m_distToXGap;
        float m_xWidth;
        float m_energy;
    };

    /**
     *  @brief Calculate all properties of a hit relevant for constructing the hit feature vector.
     *
     *  @param pCaloHit the input hit
     *  @param vtxPos the position of the 2D neutrino vertex associated with the hit
     *  @param hitFeatures struct to store calculated hit properties
     */
    static void CalculateHitFeatures(const pandora::CaloHit *const pCaloHit, const std::set<float> &detXGaps,
        pandora::CartesianVector vtxPos, HitFeatures &hitFeatures);
};

} // namespace lar_dl_content

#endif // #ifndef LAR_DL_SHOWER_HELPER_H
