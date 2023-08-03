/**
 *  @file   larpandoracontent/LArMonitoring/VertexDrawingAlgorithm.h
 *
 *  @brief  Header file for the vertex drawing algorithm.h
 *
 *  $Log: $
 */
#ifndef LAR_VERTEX_DRAWING_ALGORITHM_H
#define LAR_VERTEX_DRAWING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

/**
 *  @brief  VertexDrawingAlgorithm class
 */
class VertexDrawingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    VertexDrawingAlgorithm();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_nuVertexListName;                       ///< Name of input neutrino vertex list
};

} // namespace lar_content

#endif // LAR_VERTEX_DRAWING_ALGORITHM_H
