/**
 *  @file   larpandoracontent/LArHelpers/LArHitWidthHelper.cc
 *
 *  @brief  Implementation of the lar hit width helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar_content
{

void LArHitWidthHelper::GetExtremalCoordinatesX(const Cluster *const pCluster, CartesianVector &lowerXCoordinate, CartesianVector &higherXCoordinate, const float maxConstituentHitWidth)
{

  CartesianPointVector coordinateVector(GetConstituentHits(pCluster, maxConstituentHitWidth));

    // Find the extremal values of the X, Y and Z coordinates
    float xMin(+std::numeric_limits<float>::max());
    float yMin(+std::numeric_limits<float>::max());
    float zMin(+std::numeric_limits<float>::max());
    float xMax(-std::numeric_limits<float>::max());
    float yMax(-std::numeric_limits<float>::max());
    float zMax(-std::numeric_limits<float>::max());

    for (CartesianPointVector::const_iterator pIter = coordinateVector.begin(), pIterEnd = coordinateVector.end(); pIter != pIterEnd; ++pIter)
    {
        const CartesianVector &pos = *pIter;
        xMin = std::min(pos.GetX(), xMin);
        xMax = std::max(pos.GetX(), xMax);
        yMin = std::min(pos.GetY(), yMin);
        yMax = std::max(pos.GetY(), yMax);
        zMin = std::min(pos.GetZ(), zMin);
        zMax = std::max(pos.GetZ(), zMax);
    }

    // Choose the coordinate with the greatest span (keeping any ties)
    const float xAve(0.5f * (xMin + xMax));
    const float yAve(0.5f * (yMin + yMax));
    const float zAve(0.5f * (zMin + zMax));

    const float xSpan(xMax - xMin);
    const float ySpan(yMax - yMin);
    const float zSpan(zMax - zMin);

    const bool useX((xSpan > std::numeric_limits<float>::epsilon()) && (xSpan + std::numeric_limits<float>::epsilon() > std::max(ySpan, zSpan)));
    const bool useY((ySpan > std::numeric_limits<float>::epsilon()) && (ySpan + std::numeric_limits<float>::epsilon() > std::max(zSpan, xSpan)));
    const bool useZ((zSpan > std::numeric_limits<float>::epsilon()) && (zSpan + std::numeric_limits<float>::epsilon() > std::max(xSpan, ySpan)));

    // Find the extremal hits separately for the chosen coordinates
    CartesianPointVector candidateVector;

    for (CartesianPointVector::const_iterator pIter = coordinateVector.begin(), pIterEnd = coordinateVector.end(); pIter != pIterEnd; ++pIter)
    {
        const CartesianVector &pos = *pIter;

        if (useX)
        {
            if (((pos.GetX() - xMin) < std::numeric_limits<float>::epsilon()) || ((pos.GetX() - xMax) > -std::numeric_limits<float>::epsilon()))
                candidateVector.push_back(pos);
        }

        if (useY)
        {
            if (((pos.GetY() - yMin) < std::numeric_limits<float>::epsilon()) || ((pos.GetY() - yMax) > -std::numeric_limits<float>::epsilon()))
                candidateVector.push_back(pos);
        }

        if (useZ)
        {
            if (((pos.GetZ() - zMin) < std::numeric_limits<float>::epsilon()) || ((pos.GetZ() - zMax) > -std::numeric_limits<float>::epsilon()))
                candidateVector.push_back(pos);
        }
    }

    // Finally, find the pair of hits that are separated by the greatest distance
    CartesianVector firstCoordinate(xAve, yAve, zAve);
    CartesianVector secondCoordinate(xAve, yAve, zAve);
    float maxDistanceSquared(+std::numeric_limits<float>::epsilon());

    for (CartesianPointVector::const_iterator iterI = candidateVector.begin(), iterEndI = candidateVector.end(); iterI != iterEndI; ++iterI)
    {
        const CartesianVector &posI = *iterI;

        for (CartesianPointVector::const_iterator iterJ = iterI, iterEndJ = candidateVector.end(); iterJ != iterEndJ; ++iterJ)
        {
            const CartesianVector &posJ = *iterJ;

            const float distanceSquared((posI - posJ).GetMagnitudeSquared());

            if (distanceSquared > maxDistanceSquared)
            {
                maxDistanceSquared = distanceSquared;
                firstCoordinate = posI;
                secondCoordinate = posJ;
            }
        }
    }

    // Set the inner and outer coordinates (Check Z first, then X in the event of a tie)

    const float deltaX(secondCoordinate.GetX() - firstCoordinate.GetX());
    const float deltaZ(secondCoordinate.GetZ() - firstCoordinate.GetZ());

    if ((deltaX > 0.f) || ((std::fabs(deltaX) < std::numeric_limits<float>::epsilon()) && (deltaZ > 0.f)))
    {
        lowerXCoordinate = firstCoordinate;
        higherXCoordinate = secondCoordinate;
    }
    else
    {
        lowerXCoordinate = secondCoordinate;
        higherXCoordinate = firstCoordinate;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianPointVector LArHitWidthHelper::GetConstituentHits(const Cluster *const pCluster, const float maxConstituentHitWidth) 
{

    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    if (orderedCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    CartesianPointVector constituentHitPositionVector;

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        for(CaloHitList::const_iterator hitIter = iter->second->begin(); hitIter != iter->second->end(); ++hitIter) 
        {
            const CaloHit *const pCaloHit = *hitIter;

	    const float hitWidth = pCaloHit->GetCellSize1();
            const unsigned int numberOfConstituentHits = floor(hitWidth/maxConstituentHitWidth) + 1;
            const float constituentHitWidth = hitWidth/static_cast<float>(numberOfConstituentHits);
	 
	    // start at end of cluster with the lowest x value
            float xPositionAlongHit(pCaloHit->GetPositionVector().GetX() - (hitWidth/2));
	    for(unsigned int i(0); i < numberOfConstituentHits; ++i) 
	    {
                i == 0 ? xPositionAlongHit += constituentHitWidth/2 : xPositionAlongHit += constituentHitWidth;
	        CartesianVector consituentHitPosition(xPositionAlongHit, 0, pCaloHit->GetPositionVector().GetZ());
                constituentHitPositionVector.push_back(consituentHitPosition);
	    }
	}
    }

    std::sort(constituentHitPositionVector.begin(), constituentHitPositionVector.end(), LArClusterHelper::SortCoordinatesByPosition);
    return constituentHitPositionVector;

}

//------------------------------------------------------------------------------------------------------------------------------------------


CartesianPointVector LArHitWidthHelper::GetUniformConstituentHits(const Cluster *const pCluster, const float constituentHitWidth)
{

    OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    //std::cout << "NUMBER OF HITS: " << orderedCaloHitList.size() << std::endl;

    if (orderedCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    CartesianPointVector constituentHitPositionVector;

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        for(CaloHitList::const_iterator hitIter = iter->second->begin(); hitIter != iter->second->end(); ++hitIter) 
        {
            const CaloHit *const pCaloHit = *hitIter;

            const unsigned int numberOfConstituents(floor(pCaloHit->GetCellSize1()/constituentHitWidth) + 1);

	    //std::cout << "HIT WIDTH: " << pCaloHit->GetCellSize1() << std::endl;
	    //std::cout << "HITS TO BE BROKEN INTO: " << numberOfConstituents << std::endl;
	    //std::cout << "HIT CENTRE: " << pCaloHit->GetPositionVector() << std::endl;

            CartesianVector positionAlongHit(pCaloHit->GetPositionVector());

            if(numberOfConstituents % 2 == 1)
            {
	        constituentHitPositionVector.push_back(positionAlongHit);
	        for(unsigned int i(0); i < ((numberOfConstituents - 1)/2); ++i)
	        {
	            positionAlongHit += CartesianVector(constituentHitWidth,0,0);
		    //std::cout << "CONSTITUENT HIT PLUS ODD: " << positionAlongHit << std::endl;
		    constituentHitPositionVector.push_back(positionAlongHit);
		}

	        positionAlongHit = pCaloHit->GetPositionVector();
	        for(unsigned int i(0); i < ((numberOfConstituents - 1)/2); ++i)
	        {
	            positionAlongHit -= CartesianVector(constituentHitWidth,0,0);
                    //std::cout << "CONSTITUENT HIT MINUS ODD: " << positionAlongHit << std::endl;
		    constituentHitPositionVector.push_back(positionAlongHit);
		}
	    }
            else
	    {
                for(unsigned int i(0); i < (numberOfConstituents/2); ++i)
	        {
	            i == 0 ? positionAlongHit += CartesianVector(constituentHitWidth/2,0,0) : CartesianVector(constituentHitWidth,0,0);                
		    //std::cout << "CONSTITUENT HIT PLUS EVEN: " << positionAlongHit << std::endl;
		    constituentHitPositionVector.push_back(positionAlongHit);
		}

		positionAlongHit = pCaloHit->GetPositionVector();
	        for(unsigned int i(0); i < (numberOfConstituents/2); ++i)
	        {
	            i == 0 ? positionAlongHit -= CartesianVector(constituentHitWidth/2,0,0) : CartesianVector(constituentHitWidth,0,0);
                    //std::cout << "CONSTITUENT HIT MINUS EVEN: " << positionAlongHit << std::endl;
		    constituentHitPositionVector.push_back(positionAlongHit);
		}
	    }
	}
    }


    return constituentHitPositionVector;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArHitWidthHelper::GetExtremalCoordinatesLowerX(const Cluster *const pCluster, const float maxConstituentHitWidth) 
{

  CartesianVector lowerXCoordinate(0,0,0);
  CartesianVector higherXCoordinate(0,0,0);

  LArHitWidthHelper::GetExtremalCoordinatesX(pCluster, lowerXCoordinate, higherXCoordinate, maxConstituentHitWidth);

  return lowerXCoordinate;

}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArHitWidthHelper::GetExtremalCoordinatesHigherX(const Cluster *const pCluster, const float maxConstituentHitWidth) 
{

  CartesianVector lowerXCoordinate(0,0,0);
  CartesianVector higherXCoordinate(0,0,0);

  LArHitWidthHelper::GetExtremalCoordinatesX(pCluster, lowerXCoordinate, higherXCoordinate, maxConstituentHitWidth);

  return higherXCoordinate;

}

//------------------------------------------------------------------------------------------------------------------------------------------



float LArHitWidthHelper::GetTotalClusterWeight(const Cluster *const pCluster)
{
  
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    float hitWeight(0.0);
    
    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(); hitIter != iter->second->end(); ++hitIter)
        {
	    const CaloHit *const hit = (*hitIter);
	    hitWeight += hit->GetCellSize1();
        }
    }

    return hitWeight;

}


//------------------------------------------------------------------------------------------------------------------------------------------


bool LArHitWidthHelper::SortByMaxX(const Cluster *const pLhs, const Cluster *const pRhs, const float maxConstituentHitWidth)
{

    CartesianVector lhsLowerXEdge(0,0,0);
    CartesianVector lhsUpperXEdge(0,0,0);

    CartesianVector rhsLowerXEdge(0,0,0);
    CartesianVector rhsUpperXEdge(0,0,0);

    LArHitWidthHelper::GetExtremalCoordinatesX(pLhs, lhsLowerXEdge, lhsUpperXEdge, maxConstituentHitWidth);
    LArHitWidthHelper::GetExtremalCoordinatesX(pRhs, rhsLowerXEdge, rhsUpperXEdge, maxConstituentHitWidth);
  
    return (lhsUpperXEdge.GetX() < rhsUpperXEdge.GetX());
}


//------------------------------------------------------------------------------------------------------------------------------------------



} // namespace lar_content
