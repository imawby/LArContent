/**
 *  @file   larpandoracontent/LArTwoDReco/HitWidthClusterMergingAlgorithm.cc
 *
 *  @brief  Implementation of the hit width cluster merging algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArTwoDReco/HitWidthClusterMergingAlgorithm.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{

HitWidthClusterMergingAlgorithm::HitWidthClusterMergingAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitWidthClusterMergingAlgorithm::Run()
{
  
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    PandoraMonitoringApi::Create(this->GetPandora());
    
    for(const Cluster *const pCluster : *pClusterList) 
    {
      CartesianVector innerCoordinate(0, 0, 0);
      CartesianVector outerCoordinate(0, 0, 0);
      this->GetExtremalCoordinates(pCluster, innerCoordinate, outerCoordinate);

      PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &innerCoordinate, "Inner Coordinate", RED, 2);
      PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &outerCoordinate, "Outer Coordinate", BLUE, 2);

      PandoraMonitoringApi::Pause(this->GetPandora());
    }
    
    
    return STATUS_CODE_SUCCESS;

}


void HitWidthClusterMergingAlgorithm::GetExtremalCoordinates(const Cluster *const pCluster, CartesianVector &innerCoordinate, CartesianVector &outerCoordinate) 
{
return GetExtremalCoordinates(pCluster->GetOrderedCaloHitList(), innerCoordinate, outerCoordinate);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetExtremalCoordinates(const OrderedCaloHitList &orderedCaloHitList, CartesianVector &innerCoordinate, CartesianVector &outerCoordinate) 
{
    if (orderedCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    CartesianPointVector coordinateVector;

    // OrderedCaloHitList is a map between the pseudolayer and a calo hit list (presumably hits within that pseudolayer)
    // therefore, loop over the pseudolayers and then loop over the hits in corresponding calo hit list 
    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(); hitIter != iter->second->end(); ++hitIter)
        {
            const CaloHit *const pCaloHit = *hitIter;

	    // Take into account the width of the hits (in the x direction)
	    const CartesianVector hitCentre = pCaloHit->GetPositionVector();
	    const float hitWidth = pCaloHit->GetCellSize1();

	    CartesianVector hitUpperLimit(hitCentre.GetX() + hitWidth/2, hitCentre.GetY(), hitCentre.GetZ());
	    CartesianVector hitLowerLimit(hitCentre.GetX() - hitWidth/2, hitCentre.GetY(), hitCentre.GetZ());

            coordinateVector.push_back(hitUpperLimit);
	    coordinateVector.push_back(hitLowerLimit);
        }
    }
    std::sort(coordinateVector.begin(), coordinateVector.end(), LArClusterHelper::SortCoordinatesByPosition);
    return LArClusterHelper::GetExtremalCoordinates(coordinateVector, innerCoordinate, outerCoordinate);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitWidthClusterMergingAlgorithm::ReadSettings(const TiXmlHandle /*xmlHandle*/)
{

    return STATUS_CODE_SUCCESS;
}
  
} // namespace lar_content

