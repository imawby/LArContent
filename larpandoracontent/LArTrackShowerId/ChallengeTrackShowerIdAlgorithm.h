/**
 *  @file   larpandoracontent/LArCheating/ChallengeTrackShowerIdAlgorithm.h
 *
 *  @brief  Header file for the challenge track shower if class.
 *
 *  $Log: $
 */
#ifndef LAR_CHALLENGE_TRACK_SHOWER_ID_ALGORITHM_H
#define LAR_CHALLENGE_TRACK_SHOWER_ID_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include <unordered_map>

namespace lar_content
{

/**
 *  @brief  ChallengeTrackShowerIdAlgorithm class
 */
class ChallengeTrackShowerIdAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    ChallengeTrackShowerIdAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

};

} // namespace lar_content

#endif // #ifndef LAR_CHALLENGE_TRACK_SHOWER_ID_ALGORITHM_H
