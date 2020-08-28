/**
 *  @file   larpandoracontent/LArMonitoring/MCDeltaRayAlgorithm.h
 *
 *  @brief  Header file for the mc delta ray algorithm.
 *
 *  $Log: $
 */
#ifndef LAR_MC_DELTA_RAY_ALGORITHM_H
#define LAR_MC_DELTA_RAY_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

namespace lar_content
{

/**
 *  @brief  MCDeltaRayAlgorithm class
 */
class MCDeltaRayAlgorithm: public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    MCDeltaRayAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;
    std::string m_mcParticleListName;
};

} // namespace lar_content

#endif // LAR_MC_DELTA_RAY_ALGORITHM_H
