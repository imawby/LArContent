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
#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"
#include "Helpers/ClusterFitHelper.h"


using namespace pandora;

namespace lar_content
{


HitWidthClusterMergingAlgorithm::HitWidthClusterMergingAlgorithm() :
  m_clusterListName(),
  m_minClusterWeight(0.5),      
  m_maxXMergeDistance(5),     
  m_maxZMergeDistance(2),     
  m_maxMergeCosOpeningAngle(0.97), 
  m_maxConstituentHitWidth(0.5),
  m_fitToFullCluster(true),
  m_fittingSampleWeight(5),
  m_clusterToFitParametersMap()
{
}


//------------------------------------------------------------------------------------------------------------------------------------------

HitWidthClusterMergingAlgorithm::ClusterFitParameters::ClusterFitParameters(const pandora::Cluster *const pCluster, const float maxConstituentHitWidth) :
    m_pCluster(pCluster),
    m_numCaloHits(pCluster->GetNCaloHits()), 
    m_totalWeight(LArHitWidthHelper::GetTotalClusterWeight(pCluster)),
    m_lowerXExtrema(LArHitWidthHelper::GetExtremalCoordinatesLowerX(pCluster, maxConstituentHitWidth)), //sadly this means that the cluster break up has to be repeated
    m_higherXExtrema(LArHitWidthHelper::GetExtremalCoordinatesHigherX(pCluster, maxConstituentHitWidth)), //sadly this means that the cluster break up has to be repeated 
    m_currentClusterSort(m_higherXExtrema),
    m_testClusterSort(m_lowerXExtrema),
    m_currentClusterPositionToWeightMap(m_currentClusterSort),
    m_testClusterPositionToWeightMap(m_testClusterSort)
{

    OrderedCaloHitList orderedCaloHitList = pCluster->GetOrderedCaloHitList();

    for(OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter !=  orderedCaloHitList.end(); ++iter) 
    {
        for(CaloHitList::const_iterator hitIter = iter->second->begin(); hitIter != iter->second->end(); ++hitIter) 
        {
	    const CaloHit *const hit = (*hitIter);

	    float hitWidth = hit->GetCellSize1();
	    float hitWeight = hitWidth;

            unsigned int numberOfConstituentHits = floor(hitWidth/maxConstituentHitWidth) + 1;
            float constituentHitWidth = hitWidth/numberOfConstituentHits;
	    float constituentHitWeight = hitWeight/numberOfConstituentHits;

	    float xPositionAlongHit(hit->GetPositionVector().GetX() - (hitWidth/2));
	    for(unsigned int i(0); i < numberOfConstituentHits; ++i) 
	    {
                i == 0 ? xPositionAlongHit += constituentHitWidth/2 : xPositionAlongHit += constituentHitWidth;
	         m_currentClusterPositionToWeightMap.insert(std::pair(CartesianVector(xPositionAlongHit, 0, hit->GetPositionVector().GetZ()), constituentHitWeight));
	         m_testClusterPositionToWeightMap.insert(std::pair(CartesianVector(xPositionAlongHit, 0, hit->GetPositionVector().GetZ()), constituentHitWeight));
	    }
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{

    for (const Cluster *const pCluster : *pClusterList)
    {
        if(LArHitWidthHelper::GetTotalClusterWeight(pCluster) < m_minClusterWeight)
            continue;

        clusterVector.push_back(pCluster);
    }


    //ORDER BY MAX EXTREMAL X COORDINATE
    //HAVE TO REPEAT BREAKING UP CLUSTER IN THIS HORRIBLE WAY BECAUSE GETLISTOFCLEANCLUSTER FUNCTION IS CONSTANT
    std::sort(clusterVector.begin(), clusterVector.end(), LArHitWidthHelper::MaxXPositionSort(m_maxConstituentHitWidth));
}

//------------------------------------------------------------------------------------------------------------------------------------------


void HitWidthClusterMergingAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{

    ClusterToFitParametersMap clusterToFitParametersMap;

    // Fill cluster to fit parameter map
    for(const Cluster *const pCluster : clusterVector) 
    {
      clusterToFitParametersMap.insert(std::pair(pCluster, ClusterFitParameters(pCluster, m_maxConstituentHitWidth)));
    }
    //PandoraMonitoringApi::Create(this->GetPandora());
    //PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);

    // ATTN This method assumes that clusters have been sorted by x position (low x -> high x)
    for (ClusterVector::const_iterator iterCurrentCluster = clusterVector.begin(); iterCurrentCluster != clusterVector.end(); ++iterCurrentCluster)
    {

        const Cluster *const pCurrentCluster = *iterCurrentCluster;
	const ClusterFitParameters currentFitParameters = clusterToFitParametersMap.at(pCurrentCluster);
	
        //ClusterList currentClusterList;
        //currentClusterList.push_back(currentFitParameters.m_pCluster);
        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentClusterList, "CURRENT", BLACK, 2);



        for (ClusterVector::const_iterator iterTestCluster = iterCurrentCluster; iterTestCluster != clusterVector.end(); ++iterTestCluster)
        {
	    
            if (iterCurrentCluster == iterTestCluster)
                continue;

	    const Cluster *const pTestCluster = *iterTestCluster;
            const ClusterFitParameters testFitParameters = clusterToFitParametersMap.at(pTestCluster);


            if (!this->AreClustersAssociated(currentFitParameters, testFitParameters))
	        continue;

            clusterAssociationMap[pCurrentCluster].m_forwardAssociations.insert(pTestCluster);
            clusterAssociationMap[pTestCluster].m_backwardAssociations.insert(pCurrentCluster);
        }

        //PandoraMonitoringApi::ViewEvent(this->GetPandora());

    }

    // remove all 'shortcut' routes
    this->CleanupClusterAssociations(clusterVector, clusterAssociationMap);

}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::CleanupClusterAssociations(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{

  // Create temporary map so can delete elements whilst still iterating over them
  ClusterAssociationMap tempMap(clusterAssociationMap);

  for(const Cluster *const pCluster : clusterVector) 
  {
      const ClusterAssociationMap::const_iterator primaryMapIter = clusterAssociationMap.find(pCluster);

      if(primaryMapIter == clusterAssociationMap.end())
          continue;

      const ClusterSet &primaryForwardAssociations(primaryMapIter->second.m_forwardAssociations);

      // loop through the primary associations of each cluster
      // remove clusters that are present in secondary associations of other primary clusters 
      for(const Cluster *const pConsideredCluster : primaryForwardAssociations)
      {
	  for(const Cluster *const pPrimaryCluster : primaryForwardAssociations)
          {
	      if(pConsideredCluster == pPrimaryCluster)
	          continue;

              const ClusterAssociationMap::const_iterator secondaryMapIter = clusterAssociationMap.find(pPrimaryCluster);

	      // if primary cluster has no associations (this shouldn't ever be the case)
              if(secondaryMapIter == clusterAssociationMap.end()) 
                  continue;

	      const ClusterSet &secondaryForwardAssociations(secondaryMapIter->second.m_forwardAssociations);
	      
	      // if cluster is present in the forward associations of any cluster at the same level remove
	      if(secondaryForwardAssociations.find(pConsideredCluster) != secondaryForwardAssociations.end()) 
              {
		  ClusterSet &tempPrimaryForwardAssociations(tempMap.find(pCluster)->second.m_forwardAssociations);
		  const ClusterSet::const_iterator forwardAssociationToRemove(tempPrimaryForwardAssociations.find(pConsideredCluster));

		  // if association has already been removed
		  if(forwardAssociationToRemove == tempPrimaryForwardAssociations.end())
		      continue;

		  ClusterSet &tempPrimaryBackwardAssociations(tempMap.find(pConsideredCluster)->second.m_backwardAssociations);
		  const ClusterSet::const_iterator backwardAssociationToRemove(tempPrimaryBackwardAssociations.find(pCluster));

		  // if association has already been removed
		  if(backwardAssociationToRemove == tempPrimaryBackwardAssociations.end())
		      continue;

                  tempPrimaryForwardAssociations.erase(forwardAssociationToRemove);
		  tempPrimaryBackwardAssociations.erase(backwardAssociationToRemove);
      	      }
	      
	  }
      }
  }
  /*
  PandoraMonitoringApi::Create(this->GetPandora());
  PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);

  for(const Cluster *const pCluster : clusterVector) 
  {
      const ClusterAssociationMap::const_iterator oldIter(clusterAssociationMap.find(pCluster));
      const ClusterAssociationMap::const_iterator newIter(tempMap.find(pCluster));

      if(oldIter == clusterAssociationMap.end() || newIter == clusterAssociationMap.end())
	continue;

      ClusterList currentCluster;
      currentCluster.push_back(pCluster);
      PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "CURRENT", BLACK);

      const ClusterSet &oldForwardAssociations(oldIter->second.m_forwardAssociations);
      const ClusterSet &newForwardAssociations(newIter->second.m_forwardAssociations);

      const ClusterSet &oldBackwardAssociations(oldIter->second.m_backwardAssociations);
      const ClusterSet &newBackwardAssociations(newIter->second.m_backwardAssociations);

      ClusterList oldForwardClusters;
      for(const Cluster *const pOldCluster : oldForwardAssociations) 
      {
          oldForwardClusters.push_back(pOldCluster);
      }
      PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &oldForwardClusters, "Old Forward Clusters", DARKGREEN);

      ClusterList newForwardClusters;
      for(const Cluster *const pNewCluster : newForwardAssociations) 
      {
          newForwardClusters.push_back(pNewCluster);
      }
      PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &newForwardClusters, "New Forward Clusters", GREEN);


      ClusterList oldBackwardClusters;
      for(const Cluster *const pOldCluster : oldBackwardAssociations) 
      {
          oldBackwardClusters.push_back(pOldCluster);
      }
      PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &oldBackwardClusters, "Old Backward Clusters", DARKRED);

      ClusterList newBackwardClusters;
      for(const Cluster *const pNewCluster : newBackwardAssociations) 
      {
          newBackwardClusters.push_back(pNewCluster);
      }
      PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &newBackwardClusters, "New Backward Clusters", RED);

      PandoraMonitoringApi::ViewEvent(this->GetPandora());

  }
  */
  
  clusterAssociationMap = tempMap;

}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HitWidthClusterMergingAlgorithm::AreClustersAssociated(const ClusterFitParameters &currentFitParameters, const ClusterFitParameters &testFitParameters) const
{
  /*
    CartesianVector currentClusterDirection(GetClusterDirection(currentFitParameters.m_currentClusterPositionToWeightMap, currentFitParameters.m_numCaloHits));
    CartesianVector testClusterDirection(GetClusterDirection(testFitParameters.m_testClusterPositionToWeightMap, testFitParameters.m_numCaloHits));
    
    std::string stringTag = "TEST: " + std::to_string(std::fabs(currentClusterDirection.GetCosOpeningAngle(testClusterDirection)));


    ClusterList testClusterList;
    testClusterList.push_back(testFitParameters.m_pCluster);

    if(testFitParameters.m_lowerXExtrema.GetX() > (currentFitParameters.m_higherXExtrema.GetX() + m_maxXMergeDistance))	    
    {
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, "TEST - Outside X", RED);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testFitParameters.m_lowerXExtrema, "TEST", RED, 2);
        //PandoraMonitoringApi::Pause(this->GetPandora());
        return false;
    }

    if(testFitParameters.m_lowerXExtrema.GetZ() > (currentFitParameters.m_higherXExtrema.GetZ() + m_maxZMergeDistance) || testFitParameters.m_lowerXExtrema.GetZ() < (currentFitParameters.m_higherXExtrema.GetZ() - m_maxZMergeDistance))
    {
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, "TEST - Outside Z", RED);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testFitParameters.m_lowerXExtrema, "TEST", RED, 2);
        //PandoraMonitoringApi::Pause(this->GetPandora());
        return false;
    }

    if(fabs(currentClusterDirection.GetCosOpeningAngle(testClusterDirection)) < m_maxMergeCosOpeningAngle)
    {
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, stringTag, RED);
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testFitParameters.m_lowerXExtrema, "TEST", RED, 2);
        //PandoraMonitoringApi::Pause(this->GetPandora());
        return false;
    }
    
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testFitParameters.m_lowerXExtrema, "TEST", GREEN, 2);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, stringTag, GREEN);
    
    //PandoraMonitoringApi::Pause(this->GetPandora());

    return true;
    /*/
    

    CartesianVector currentClusterDirection(GetClusterDirection(currentFitParameters.m_currentClusterPositionToWeightMap, currentFitParameters.m_numCaloHits));
    CartesianVector testClusterDirection(GetClusterDirection(testFitParameters.m_testClusterPositionToWeightMap, testFitParameters.m_numCaloHits));

    if(testFitParameters.m_lowerXExtrema.GetX() > (currentFitParameters.m_higherXExtrema.GetX() + m_maxXMergeDistance)) {  
      return false;
    }

    if(testFitParameters.m_lowerXExtrema.GetZ() > (currentFitParameters.m_higherXExtrema.GetZ() + m_maxZMergeDistance) || testFitParameters.m_lowerXExtrema.GetZ() < (currentFitParameters.m_higherXExtrema.GetZ() - m_maxZMergeDistance)) {
      return false;
    }
    
    if(fabs(currentClusterDirection.GetCosOpeningAngle(testClusterDirection)) < m_maxMergeCosOpeningAngle) {
      return false;
    }
    
    return true;
    

}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HitWidthClusterMergingAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster,  const Cluster *const pTestCluster) const
{
  
    //NEED TO IMPLEMENT USING THE MAP INSTEAD

    CartesianVector currentMin(0,0,0), currentMax(0,0,0);
    LArHitWidthHelper::GetExtremalCoordinatesX(pCurrentCluster, currentMin, currentMax, m_maxConstituentHitWidth);

    CartesianVector testMin(0,0,0), testMax(0,0,0);
    LArHitWidthHelper::GetExtremalCoordinatesX(pTestCluster, testMin, testMax, m_maxConstituentHitWidth);

    float testMaxX(testMax.GetX()), currentMaxX(currentMax.GetX());

    if (isForward)
    {
        if (std::fabs(testMaxX - currentMaxX) > std::numeric_limits<float>::epsilon())
            return (testMaxX > currentMaxX);
    }
    else
    {
        if (std::fabs(testMaxX - currentMaxX) > std::numeric_limits<float>::epsilon())
            return (testMaxX < currentMaxX);
    }
  
    return LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster);

}


//------------------------------------------------------------------------------------------------------------------------------------------


  CartesianVector HitWidthClusterMergingAlgorithm::GetClusterDirection(const ClusterFitParameters::ClusterPositionToWeightMap &clusterPositionToWeightMap, unsigned int clusterCaloHits) const 
{

    // If cluster composed of one hit, return a transverse line as fit
    if(clusterCaloHits == 1) {
      return CartesianVector(1, 0, 0);
    }
   

    // minimise the longitudinal distance in fit (works best for transverse fits)
    CartesianVector LSTransverseClusterFitDirection(0,0,0);
    CartesianVector LSTransverseIntercept(0,0,0);
    float LSTransverseChiSquared(0);

    // minimise the transverse coordinate in fit (works best for longitudinal fits)
    CartesianVector LSLongitudinalClusterFitDirection(0,0,0);
    CartesianVector LSLongitudinalIntercept(0,0,0);
    float LSLongitudinalChiSquared(0);
    
    GetWeightedGradient(clusterPositionToWeightMap, true, LSTransverseClusterFitDirection, LSTransverseIntercept, LSTransverseChiSquared, clusterCaloHits);    
    GetWeightedGradient(clusterPositionToWeightMap, false, LSLongitudinalClusterFitDirection, LSLongitudinalIntercept, LSLongitudinalChiSquared, clusterCaloHits);

    /*
    //CODE TO DRAW EACH OF THE FITS FOR THE CLUSTER
    
    //std::cout << "Transverse Chi Squared: " << LSTransverseChiSquared << std::endl;
    //std::cout << "Longitudinal Chi Squared: " << LSLongitudinalChiSquared << std::endl;

    if(LSTransverseChiSquared < LSLongitudinalChiSquared) {
        float clusterTransverseFitGradient = (LSTransverseClusterFitDirection.GetZ()/LSTransverseClusterFitDirection.GetX());
        CartesianVector LSTransverseStartPoint(400, 0, clusterTransverseFitGradient*(400) + LSTransverseIntercept.GetZ());
        CartesianVector LSTransverseEndPoint(-400, 0, clusterTransverseFitGradient*(-400) + LSTransverseIntercept.GetZ());
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &LSTransverseStartPoint, &LSTransverseEndPoint, "Transverse", BLUE, 2, 2);
    } else {
        float clusterLongitudinalFitGradient = (LSLongitudinalClusterFitDirection.GetZ()/LSLongitudinalClusterFitDirection.GetX());
        CartesianVector LSLongitudinalStartPoint(400, 0, clusterLongitudinalFitGradient*(400) + LSLongitudinalIntercept.GetZ());
        CartesianVector LSLongitudinalEndPoint(-400, 0, clusterLongitudinalFitGradient*(-400) + LSLongitudinalIntercept.GetZ());
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &LSLongitudinalStartPoint, &LSLongitudinalEndPoint, "Longitudinal", DARKGREEN, 2, 2);
       }

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */

    // return fit with the lowest chi-squared
    if(LSTransverseChiSquared < LSLongitudinalChiSquared)
      return LSTransverseClusterFitDirection;

    return LSLongitudinalClusterFitDirection;

}
  
//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector HitWidthClusterMergingAlgorithm::GetClusterZIntercept(const ClusterFitParameters::ClusterPositionToWeightMap &clusterPositionToWeightMap, unsigned int clusterCaloHits) const
{

    // If cluster composed of one hit, return a transverse line as fit
    if(clusterCaloHits == 1) {
      return CartesianVector(0, 0, clusterPositionToWeightMap.begin()->first.GetZ());
    }      

    // minimise the longitudinal distance in fit (works best for transverse fits)
    CartesianVector LSTransverseClusterFitDirection(0,0,0);
    CartesianVector LSTransverseIntercept(0,0,0);
    float LSTransverseChiSquared(0);

    // minimise the transverse coordinate in fit (works best for longitudinal fits)
    CartesianVector LSLongitudinalClusterFitDirection(0,0,0);
    CartesianVector LSLongitudinalIntercept(0,0,0);
    float LSLongitudinalChiSquared(0);

    GetWeightedGradient(clusterPositionToWeightMap, true, LSTransverseClusterFitDirection, LSTransverseIntercept, LSTransverseChiSquared, clusterCaloHits);
    GetWeightedGradient(clusterPositionToWeightMap, false, LSLongitudinalClusterFitDirection, LSLongitudinalIntercept, LSLongitudinalChiSquared, clusterCaloHits);

    // return fit with the lowest chi-squared
    if(LSTransverseChiSquared < LSLongitudinalChiSquared)
      return LSTransverseIntercept;

    return LSLongitudinalIntercept;

}

//------------------------------------------------------------------------------------------------------------------------------------------


  void HitWidthClusterMergingAlgorithm::GetWeightedGradient(const ClusterFitParameters::ClusterPositionToWeightMap &clusterPositionToWeightMap, bool isTransverse, CartesianVector &direction, CartesianVector &intercept, float &chiSquared, float clusterCaloHits) const
{

    // cannot make a longitudinal fit to a single hit cluster
    if(!isTransverse && clusterCaloHits == 1) 
    {
      std::cout << "WARNING - CANNOT MAKE LONGITUDINAL FIT TO SINGLE HIT CLUSTER, RETURNED" << std::endl;
      throw;
    }

    float weightSum(0);
    float weightedXSum(0);
    float weightedZSum(0);

    float weightCount(0);

    bool isXConstant(true);
    bool isZConstant(true);

    // Include in fit hits nearest the relevant x extrema
    // Stop including hits when cumulative weight exceeds m_fittingSampleWeight 
    // or have already included all cluster constituent hits 
    for(ClusterFitParameters::ClusterPositionToWeightMap::const_iterator iter = clusterPositionToWeightMap.begin(); iter != clusterPositionToWeightMap.end(); ++iter) 
    {

        const CartesianVector hitPosition = iter->first;
	const float hitWeight = iter->second;
	weightCount += hitWeight;

	if(std::fabs(clusterPositionToWeightMap.begin()->first.GetX() - hitPosition.GetX()) > std::numeric_limits<float>::epsilon())
	  isXConstant = false;

	if(std::fabs(clusterPositionToWeightMap.begin()->first.GetZ() - hitPosition.GetZ()) > std::numeric_limits<float>::epsilon())
	  isZConstant = false;

	weightSum += hitWeight;
        weightedXSum += hitPosition.GetX() * hitWeight;
	weightedZSum += hitPosition.GetZ() * hitWeight;

	// think its better to go over the weight limit then under the weight limit 
        if(!m_fitToFullCluster && (weightCount > m_fittingSampleWeight))
            break;
    }

    // TO FIT A STRAIGHT LINE TO CLUSTERS WITH CONSTANT X OR Z
    if(isXConstant) 
    {
      direction = CartesianVector(0, 0, 1);
      intercept = CartesianVector(0, 0, -std::numeric_limits<float>::max());
      chiSquared = 0;
      return;
    }

    if(isZConstant) 
    {
      direction = CartesianVector(1, 0, 0);
      intercept = CartesianVector(0, 0, clusterPositionToWeightMap.begin()->first.GetZ());
      chiSquared = 0;
      return;
    }

    float weightedXMean(weightedXSum/weightSum);
    float weightedZMean(weightedZSum/weightSum); 

    float numerator(0);
    float denominator(0);
    float chi(0);

    weightCount = 0;
    for(ClusterFitParameters::ClusterPositionToWeightMap::const_iterator iter = clusterPositionToWeightMap.begin(); iter != clusterPositionToWeightMap.end(); ++iter) 
    {
        const CartesianVector hitPosition = iter->first;
	const float hitWeight = iter->second;
	weightCount += hitWeight;

	numerator += hitWeight * (hitPosition.GetX() - weightedXMean) * (hitPosition.GetZ() - weightedZMean);
	isTransverse ? denominator += hitWeight * pow(hitPosition.GetX() - weightedXMean, 2) : denominator += hitWeight * pow(hitPosition.GetZ() - weightedZMean, 2);

        if(!m_fitToFullCluster && (weightCount > m_fittingSampleWeight))
            break;
    }

    float gradient = numerator/denominator;
    isTransverse ? intercept.SetValues(0, 0, weightedZMean - gradient * weightedXMean) : intercept.SetValues(weightedXMean - gradient * weightedZMean, 0, 0);

    weightCount = 0;
    for(ClusterFitParameters::ClusterPositionToWeightMap::const_iterator iter = clusterPositionToWeightMap.begin(); iter != clusterPositionToWeightMap.end(); ++iter) 
    {
        const CartesianVector hitPosition = iter->first;
	const float hitWeight = iter->second;
	weightCount += hitWeight;

	isTransverse ? chi += hitWeight*pow(hitPosition.GetZ() - intercept.GetZ() - gradient*hitPosition.GetX(), 2) : chi += hitWeight * pow(hitPosition.GetX() - intercept.GetX() - gradient * hitPosition.GetZ(), 2);

        if(!m_fitToFullCluster && (weightCount > m_fittingSampleWeight))
            break;
    }


    // change coordinates to y=mx+c fit and normalise
    isTransverse? direction = CartesianVector(1.0, 0, gradient).GetUnitVector() : direction = CartesianVector(gradient, 0, 1.0).GetUnitVector();

    if(!isTransverse)
        intercept.SetValues(0, 0, -intercept.GetX()/gradient);


    chiSquared = chi;

    return;

}



//------------------------------------------------------------------------------------------------------------------------------------------



StatusCode HitWidthClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterWeight", m_minClusterWeight));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxXMergeDistance", m_maxXMergeDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxZMergeDistance", m_maxZMergeDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxMergeCosOpeningAngle", m_maxMergeCosOpeningAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxConstituentHitWidth", m_maxConstituentHitWidth));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FitToFullCluster", m_fitToFullCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FittingSampleWeight", m_fittingSampleWeight));


    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}
  
} // namespace lar_content



//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/*
void HitWidthClusterMergingAlgorithm::TestPopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{    

    // Fill cluster fit map 
    ClusterToFitParametersMap clusterToFitParametersMap;

    for(const Cluster *const pCluster : clusterVector) 
    {
      clusterToFitParametersMap.insert(std::pair(pCluster, ClusterFitParameters(pCluster)));
    }


    // ATTN This method assumes that clusters have been sorted by x position (low x -> high x)
    for (ClusterVector::const_iterator iterCurrentCluster = clusterVector.begin(); iterCurrentCluster != clusterVector.end(); ++iterCurrentCluster)
    {

        const Cluster *const pCurrentCluster = *iterCurrentCluster;
	const ClusterFitParameters currentFitParameters = clusterToFitParametersMap.at(pCurrentCluster);
	
        ClusterList currentClusterList;
        currentClusterList.push_back(pCurrentCluster);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentClusterList, "CURRENT", BLUE);

        CartesianVector currentLowerXEdge(currentFitParameters.m_lowerXExtrema);
        CartesianVector currentUpperXEdge(currentFitParameters.m_higherXExtrema);
            
        CartesianVector mergeXDist(currentUpperXEdge.GetX() + m_maxXMergeDistance, 0, currentUpperXEdge.GetZ() - m_maxZMergeDistance);
        CartesianVector mergeZDist1(currentUpperXEdge.GetX(), 0, currentUpperXEdge.GetZ() + m_maxZMergeDistance);
        CartesianVector mergeZDist2(currentUpperXEdge.GetX(), 0, currentUpperXEdge.GetZ() - m_maxZMergeDistance);
        CartesianVector farPoint(currentUpperXEdge.GetX() + m_maxXMergeDistance, 0, currentUpperXEdge.GetZ() + m_maxZMergeDistance);

        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &mergeZDist2, &mergeZDist1, "1", RED, 2, 2);
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &mergeZDist2, &mergeXDist, "2", RED, 2, 2);
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &farPoint, &mergeZDist1, "3", RED, 2, 2);
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &farPoint, &mergeXDist, "4", RED, 2, 2);

	PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentUpperXEdge, "Upper Coordinate", BLUE, 2);

        for (ClusterVector::const_iterator iterTestCluster = iterCurrentCluster; iterTestCluster != clusterVector.end(); ++iterTestCluster)
        {
	    
            if (iterCurrentCluster == iterTestCluster)
                continue;

	    const Cluster *const pTestCluster = *iterTestCluster;
            const ClusterFitParameters testFitParameters = clusterToFitParametersMap.at(pTestCluster);
	    
            if (!this->TestAreClustersAssociated(currentFitParameters, testFitParameters, pTestCluster))
	        continue;

            clusterAssociationMap[pCurrentCluster].m_forwardAssociations.insert(pTestCluster);
            clusterAssociationMap[pTestCluster].m_backwardAssociations.insert(pCurrentCluster);
        }

	//PandoraMonitoringApi::Pause(this->GetPandora());
        PandoraMonitoringApi::ViewEvent(this->GetPandora()); 

    }
	

}

//------------------------------------------------------------------------------------------------------------------------------------------

  bool HitWidthClusterMergingAlgorithm::TestAreClustersAssociated(const ClusterFitParameters &currentFitParameters, const ClusterFitParameters &testFitParameters, const Cluster *const pCluster) const
{

    CartesianVector currentClusterDirection(GetClusterEndDirection(currentFitParameters.m_currentClusterPositionToWeightMap, currentFitParameters.m_numCaloHits));
    CartesianVector testClusterDirection(GetClusterEndDirection(testFitParameters.m_testClusterPositionToWeightMap, testFitParameters.m_numCaloHits));
    
    std::string stringTag = "TEST: " + std::to_string(std::fabs(currentClusterDirection.GetCosOpeningAngle(testClusterDirection)));

    ClusterList testClusterList;
    testClusterList.push_back(pCluster);
    //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, stringTag, GREEN);
	    
    if(testFitParameters.m_lowerXExtrema.GetX() > (currentFitParameters.m_higherXExtrema.GetX() + m_maxXMergeDistance)) {
      PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, "TEST - Outside X", RED);
      PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testFitParameters.m_lowerXExtrema, "TEST", RED, 2);
      //PandoraMonitoringApi::Pause(this->GetPandora());
      return false;
    }

    if(testFitParameters.m_lowerXExtrema.GetZ() > (currentFitParameters.m_higherXExtrema.GetZ() + m_maxZMergeDistance) || testFitParameters.m_lowerXExtrema.GetZ() < (currentFitParameters.m_higherXExtrema.GetZ() - m_maxZMergeDistance)) {
      PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, "TEST - Outside Z", RED);
      PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testFitParameters.m_lowerXExtrema, "TEST", RED, 2);
      //PandoraMonitoringApi::Pause(this->GetPandora());
      return false;
    }
    
    if(fabs(currentClusterDirection.GetCosOpeningAngle(testClusterDirection)) < m_maxMergeCosOpeningAngle) {
      PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, stringTag, RED);
      PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testFitParameters.m_lowerXExtrema, "TEST", RED, 2);
      //PandoraMonitoringApi::Pause(this->GetPandora());
      return false;
    }
    
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testFitParameters.m_lowerXExtrema, "TEST", GREEN, 2);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, stringTag, GREEN);
    
    //PandoraMonitoringApi::Pause(this->GetPandora());

    return true;
    
}
*/
