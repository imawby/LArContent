/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLSecondaryVertexSplittingAlgorithm.h
 *
 *  @brief  Header file for the DL neutrino hierarchy algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_SECONDARY_VERTEX_SPLITTING_ALGORITHM_H
#define LAR_DL_SECONDARY_VERTEX_SPLITTING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

using namespace lar_content;

namespace lar_dl_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  DLSecondaryVertexSplittingAlgorithm class
 */
class DLSecondaryVertexSplittingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    DLSecondaryVertexSplittingAlgorithm();

private:
    pandora::StatusCode Run();    

    int GetNKinkPoints(const TwoDSlidingFitResult &slidingFit, const pandora::VertexList *const pSecVtxList);

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    bool m_visualise;
    std::string m_clusterListNameU;
    std::string m_clusterListNameV;
    std::string m_clusterListNameW;
    std::string m_secondaryVertexListName;
    std::string m_reclusteringAlgorithmName; ///< Name of daughter algorithm to use for cluster re-building
    float m_proximityForMatch;


};

} // namespace lar_dl_content

#endif // #ifndef LAR_DL_SECONDARY_VERTEX_SPLITTING_ALGORITHM_H
