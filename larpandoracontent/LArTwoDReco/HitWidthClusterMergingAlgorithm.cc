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
#include "Helpers/ClusterFitHelper.h"

using namespace pandora;

namespace lar_content
{

  HitWidthClusterMergingAlgorithm::ClusterFitParameters::ClusterFitParameters(const pandora::Cluster *const pCluster) :
    m_numCaloHits(pCluster->GetNCaloHits()), 
    m_totalWeight(GetTotalClusterWeight(pCluster)),
    m_lowerXExtrema(GetExtremalCoordinatesLowerX(pCluster)), //sadly this means that the cluster break up has to be repeated
    m_higherXExtrema(GetExtremalCoordinatesHigherX(pCluster)), //sadly this means that the cluster break up has to be repeated 
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

            unsigned int numberOfConstituentHits = floor(hitWidth/hit->GetCellSize0()) + 1;
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

HitWidthClusterMergingAlgorithm::HitWidthClusterMergingAlgorithm() :
  m_clusterListName(),
  m_minClusterWeight(0),        //ATTN - THIS HAS PROBABLY BEEN TUNED IN THE XML FILE
  m_maxXMergeDistance(0),       //ATTN - THIS HAS PROBABLY BEEN TUNED IN THE XML FILE
  m_maxZMergeDistance(0),       //ATTN - THIS HAS PROBABLY BEEN TUNED IN THE XML FILE
  m_maxMergeCosOpeningAngle(0), //ATTN - THIS HAS PROBABLY BEEN TUNED IN THE XML FILE
  m_fittingSampleWeight(100),       //ATTN - THIS HAS PROBABLY BEEN TUNED IN THE XML FILE
  m_clusterToFitParametersMap()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------


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

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    //ignore clusters that aren't transverse??
    CartesianVector origin(0,0,0);
    ClusterFitParameters::ClusterPositionSort distanceToOriginSort(origin);

    for (const Cluster *const pCluster : *pClusterList)
    {
        if(GetTotalClusterWeight(pCluster) < m_minClusterWeight)
            continue;

        clusterVector.push_back(pCluster);
    }

    //ORDER BY MAX EXTREMAL X COORDINATE
    std::sort(clusterVector.begin(), clusterVector.end(), SortByMaxX);
}

//------------------------------------------------------------------------------------------------------------------------------------------


void HitWidthClusterMergingAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
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
    }

    // remove all 'shortcut' routes
    this->CleanupForwardAssociations(clusterVector, clusterAssociationMap);

}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::CleanupForwardAssociations(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{

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

    CartesianVector currentClusterDirection(GetClusterEndDirection(currentFitParameters.m_currentClusterPositionToWeightMap, currentFitParameters.m_numCaloHits));
    CartesianVector testClusterDirection(GetClusterEndDirection(testFitParameters.m_testClusterPositionToWeightMap, testFitParameters.m_numCaloHits));

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
    HitWidthClusterMergingAlgorithm::GetExtremalCoordinates(pCurrentCluster, currentMin, currentMax);

    CartesianVector testMin(0,0,0), testMax(0,0,0);
    HitWidthClusterMergingAlgorithm::GetExtremalCoordinates(pTestCluster, testMin, testMax);

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


  CartesianVector HitWidthClusterMergingAlgorithm::GetClusterEndDirection(const ClusterFitParameters::ClusterPositionToWeightMap &clusterPositionToWeightMap, unsigned int clusterCaloHits) const 
{

    if(clusterCaloHits == 1) {
      //std::cout << "Has single hit" << std::endl;
      return CartesianVector(1, 0, 0);
    }
   
    CartesianVector LSTransverseClusterFitDirection(0,0,0);
    CartesianVector LSTransverseIntercept(0,0,0);
    float LSTransverseChiSquared(0);

    CartesianVector LSLongitudinalClusterFitDirection(0,0,0);
    CartesianVector LSLongitudinalIntercept(0,0,0);
    float LSLongitudinalChiSquared(0);
    
    GetWeightedSubGradient(clusterPositionToWeightMap, true, LSTransverseClusterFitDirection, LSTransverseIntercept, LSTransverseChiSquared, clusterCaloHits);    
    GetWeightedSubGradient(clusterPositionToWeightMap, false, LSLongitudinalClusterFitDirection, LSLongitudinalIntercept, LSLongitudinalChiSquared, clusterCaloHits);

    
    //CODE TO DRAW EACH OF THE FITS FOR THE CLUSTER
    /*
    //std::cout << "Transverse Chi Squared: " << LSTransverseChiSquared << std::endl;
    //std::cout << "Longitudinal Chi Squared: " << LSLongitudinalChiSquared << std::endl;

    if(LSTransverseChiSquared < LSLongitudinalChiSquared) {
        float clusterTransverseFitGradient = (LSTransverseClusterFitDirection.GetZ()/LSTransverseClusterFitDirection.GetX());
        CartesianVector LSTransverseStartPoint(LSTransverseIntercept);
        CartesianVector LSTransverseEndPoint(-400, 0, clusterTransverseFitGradient*(-400) + LSTransverseIntercept.GetZ());
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &LSTransverseStartPoint, &LSTransverseEndPoint, "Transverse", BLUE, 2, 2);
    } else {
        float clusterLongitudinalFitGradient = (LSLongitudinalClusterFitDirection.GetZ()/LSLongitudinalClusterFitDirection.GetX());
        CartesianVector LSLongitudinalStartPoint(LSLongitudinalIntercept);
        CartesianVector LSLongitudinalEndPoint(-400, 0, clusterLongitudinalFitGradient*(-400) + LSLongitudinalIntercept.GetZ());
        PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &LSLongitudinalStartPoint, &LSLongitudinalEndPoint, "Longitudinal", DARKGREEN, 2, 2);
    }

    //PandoraMonitoringApi::Pause(this->GetPandora());
    */

    if(LSTransverseChiSquared < LSLongitudinalChiSquared)
      return LSTransverseClusterFitDirection;

    return LSLongitudinalClusterFitDirection;

}
  
//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector HitWidthClusterMergingAlgorithm::GetClusterEndZIntercept(const ClusterFitParameters::ClusterPositionToWeightMap &clusterPositionToWeightMap, unsigned int clusterCaloHits) const
{

    if(clusterCaloHits == 1) {
      std::cout << "HAS SINGLE HIT" << std::endl;
      return CartesianVector(0, 0, clusterPositionToWeightMap.begin()->first.GetZ());
    }      

    CartesianVector LSTransverseClusterFitDirection(0,0,0);
    CartesianVector LSTransverseIntercept(0,0,0);
    float LSTransverseChiSquared(0);

    CartesianVector LSLongitudinalClusterFitDirection(0,0,0);
    CartesianVector LSLongitudinalIntercept(0,0,0);
    float LSLongitudinalChiSquared(0);

    GetWeightedSubGradient(clusterPositionToWeightMap, true, LSTransverseClusterFitDirection, LSTransverseIntercept, LSTransverseChiSquared, clusterCaloHits);
    GetWeightedSubGradient(clusterPositionToWeightMap, false, LSLongitudinalClusterFitDirection, LSLongitudinalIntercept, LSLongitudinalChiSquared, clusterCaloHits);

    if(LSTransverseChiSquared < LSLongitudinalChiSquared)
      return LSTransverseIntercept;

    return LSLongitudinalIntercept;

}

//------------------------------------------------------------------------------------------------------------------------------------------


  void HitWidthClusterMergingAlgorithm::GetWeightedSubGradient(const ClusterFitParameters::ClusterPositionToWeightMap &clusterPositionToWeightMap, bool isTransverse, CartesianVector &direction, CartesianVector &intercept, float &chiSquared, float clusterCaloHits) const
{

  //std::cout << "Number of calo hits: " << clusterCaloHits << std::endl;

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

	//std::cout << "Hit position: " << hitPosition << std::endl;
        //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "Fit Point", BLUE, 2);

	if(std::fabs(clusterPositionToWeightMap.begin()->first.GetX() - hitPosition.GetX()) > std::numeric_limits<float>::epsilon())
	  isXConstant = false;

	if(std::fabs(clusterPositionToWeightMap.begin()->first.GetZ() - hitPosition.GetZ()) > std::numeric_limits<float>::epsilon())
	  isZConstant = false;

	weightSum += hitWeight;
        weightedXSum += hitPosition.GetX() * hitWeight;
	weightedZSum += hitPosition.GetZ() * hitWeight;

	// think its better to go over the weight limit then under the weight limit 
        if(weightCount > m_fittingSampleWeight)
            break;
    }

    //std::cout << "Weight included in fit: " << weightSum << std::endl;

    // TO FIT A STRAIGHT LINE TO CLUSTERS WITH CONSTANT X OR Z
    if(isXConstant) 
    {
      //std::cout << "CONSTANT X" << std::endl;
      direction = CartesianVector(0, 0, 1);
      intercept = CartesianVector(0, 0, -std::numeric_limits<float>::max());
      chiSquared = 0;
      return;
    }

    if(isZConstant) 
    {
      //std::cout << "CONSTANT Z" << std::endl;
      direction = CartesianVector(1, 0, 0);
      intercept = CartesianVector(0, 0, clusterPositionToWeightMap.begin()->first.GetZ());
      chiSquared = 0;
      return;
    }


    //std::cout << "WeightedXSum: " << weightedXSum << std::endl;
    //std::cout << "WeightedZSum: " << weightedZSum << std::endl;
    //std::cout << "WeightSum: " << weightSum << std::endl;

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

        if(weightCount > m_fittingSampleWeight)
            break;
 
    }

    //std::cout << "Denominator: " << denominator << std::endl;
    //std::cout << "Numerator: " << numerator << std::endl;

    float gradient = numerator/denominator;
    isTransverse ? intercept.SetValues(0, 0, weightedZMean - gradient * weightedXMean) : intercept.SetValues(weightedXMean - gradient * weightedZMean, 0, 0);

    //std::cout << "Gradient: " << gradient << std::endl;
    //std::cout << "Intercept: " << intercept << std::endl;

    weightCount = 0;
    for(ClusterFitParameters::ClusterPositionToWeightMap::const_iterator iter = clusterPositionToWeightMap.begin(); iter != clusterPositionToWeightMap.end(); ++iter) 
    {

        const CartesianVector hitPosition = iter->first;
	const float hitWeight = iter->second;
	weightCount += hitWeight;


	isTransverse ? chi += hitWeight*pow(hitPosition.GetZ() - intercept.GetZ() - gradient*hitPosition.GetX(), 2) : chi += hitWeight * pow(hitPosition.GetX() - intercept.GetX() - gradient * hitPosition.GetZ(), 2);

        if(weightCount > m_fittingSampleWeight)
            break;

    }


    isTransverse? direction = CartesianVector(1.0, 0, gradient).GetUnitVector() : direction = CartesianVector(gradient, 0, 1.0).GetUnitVector();


    if(!isTransverse)
        intercept.SetValues(0, 0, -intercept.GetX()/gradient);


    chiSquared = chi;

    return;

}

//------------------------------------------------------------------------------------------------------------------------------------------

  CartesianVector HitWidthClusterMergingAlgorithm::GetClusterDirection(const ClusterFitParameters::ClusterPositionToWeightMap &clusterPositionToWeightMap, unsigned int clusterCaloHits) const 
{

    // If has single hit, put cluster direction in the direction of x
    if(clusterCaloHits == 1) {
      return CartesianVector(1, 0, 0);
    }
   
    CartesianVector LSTransverseClusterFitDirection(0,0,0);
    CartesianVector LSTransverseIntercept(0,0,0);
    float LSTransverseChiSquared(0);

    CartesianVector LSLongitudinalClusterFitDirection(0,0,0);
    CartesianVector LSLongitudinalIntercept(0,0,0);
    float LSLongitudinalChiSquared(0);
    
    GetWeightedGradient(clusterPositionToWeightMap, true, LSTransverseClusterFitDirection, LSTransverseIntercept, LSTransverseChiSquared, clusterCaloHits);    
    GetWeightedGradient(clusterPositionToWeightMap, false, LSLongitudinalClusterFitDirection, LSLongitudinalIntercept, LSLongitudinalChiSquared, clusterCaloHits);
    
    if(LSTransverseChiSquared < LSLongitudinalChiSquared)
      return LSTransverseClusterFitDirection;

    return LSLongitudinalClusterFitDirection;

}
  
//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector HitWidthClusterMergingAlgorithm::GetClusterZIntercept(const ClusterFitParameters::ClusterPositionToWeightMap &clusterPositionToWeightMap, unsigned int clusterCaloHits) const
{

    if(clusterCaloHits == 1) {
      return CartesianVector(0, 0, clusterPositionToWeightMap.begin()->first.GetZ());
    }      

    CartesianVector LSTransverseClusterFitDirection(0,0,0);
    CartesianVector LSTransverseIntercept(0,0,0);
    float LSTransverseChiSquared(0);

    CartesianVector LSLongitudinalClusterFitDirection(0,0,0);
    CartesianVector LSLongitudinalIntercept(0,0,0);
    float LSLongitudinalChiSquared(0);

    GetWeightedGradient(clusterPositionToWeightMap, true, LSTransverseClusterFitDirection, LSTransverseIntercept, LSTransverseChiSquared, clusterCaloHits);
    GetWeightedGradient(clusterPositionToWeightMap, false, LSLongitudinalClusterFitDirection, LSLongitudinalIntercept, LSLongitudinalChiSquared, clusterCaloHits);

    if(LSTransverseChiSquared < LSLongitudinalChiSquared)
      return LSTransverseIntercept;

    return LSLongitudinalIntercept;

}

//------------------------------------------------------------------------------------------------------------------------------------------


  void HitWidthClusterMergingAlgorithm::GetWeightedGradient(const ClusterFitParameters::ClusterPositionToWeightMap &clusterPositionToWeightMap, bool isTransverse, CartesianVector &direction, CartesianVector &intercept, float &chiSquared, float clusterCaloHits) const
{

    if(!isTransverse && clusterCaloHits == 1) 
    {
      std::cout << "WARNING - CANNOT MAKE LONGITUDINAL FIT TO SINGLE HIT CLUSTER, RETURNED" << std::endl;
      throw;
    }

    float weightSum(0);
    float weightedXSum(0);
    float weightedZSum(0);

    bool isXConstant(true);
    bool isZConstant(true);

    for(ClusterFitParameters::ClusterPositionToWeightMap::const_iterator iter = clusterPositionToWeightMap.begin(); iter != clusterPositionToWeightMap.end(); ++iter) 
    {

        const CartesianVector hitPosition = iter->first;
	const float hitWeight = iter->second;

	if(std::fabs(clusterPositionToWeightMap.begin()->first.GetX() - hitPosition.GetX()) > std::numeric_limits<float>::epsilon())
	  isXConstant = false;

	if(std::fabs(clusterPositionToWeightMap.begin()->first.GetZ() - hitPosition.GetZ()) > std::numeric_limits<float>::epsilon())
	  isZConstant = false;

	weightSum += hitWeight;
        weightedXSum += hitPosition.GetX() * hitWeight;
	weightedZSum += hitPosition.GetZ() * hitWeight;

    }

    // TO FIT A STRAIGHT LINE TO CLUSTERS WITH CONSTANT X OR Z
    if(isXConstant) 
    {
      //std::cout << "CONSTANT X" << std::endl;
      direction = CartesianVector(0, 0, 1);
      intercept = CartesianVector(0, 0, -std::numeric_limits<float>::max());
      chiSquared = 0;
      return;
    }

    if(isZConstant) 
    {
      //std::cout << "CONSTANT Z" << std::endl;
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

    for(ClusterFitParameters::ClusterPositionToWeightMap::const_iterator iter = clusterPositionToWeightMap.begin(); iter != clusterPositionToWeightMap.end(); ++iter) 
    {
        const CartesianVector hitPosition = iter->first;
	const float hitWeight = iter->second;

	numerator += hitWeight * (hitPosition.GetX() - weightedXMean) * (hitPosition.GetZ() - weightedZMean);
	isTransverse ? denominator += hitWeight * pow(hitPosition.GetX() - weightedXMean, 2) : denominator += hitWeight * pow(hitPosition.GetZ() - weightedZMean, 2);

    }

    float gradient = numerator/denominator;
    isTransverse ? intercept.SetValues(0, 0, weightedZMean - gradient * weightedXMean) : intercept.SetValues(weightedXMean - gradient * weightedZMean, 0, 0);


    for(ClusterFitParameters::ClusterPositionToWeightMap::const_iterator iter = clusterPositionToWeightMap.begin(); iter != clusterPositionToWeightMap.end(); ++iter) 
    {

        const CartesianVector hitPosition = iter->first;
	const float hitWeight = iter->second;


	isTransverse ? chi += hitWeight*pow(hitPosition.GetZ() - intercept.GetZ() - gradient*hitPosition.GetX(), 2) : chi += hitWeight * pow(hitPosition.GetX() - intercept.GetX() - gradient * hitPosition.GetZ(), 2);

    }

    isTransverse? direction = CartesianVector(1.0, 0, gradient).GetUnitVector() : direction = CartesianVector(gradient, 0, 1.0).GetUnitVector();


    if(!isTransverse)
        intercept.SetValues(0, 0, -intercept.GetX()/gradient);

    chiSquared = chi;

    return;

}


//------------------------------------------------------------------------------------------------------------------------------------------


CartesianVector HitWidthClusterMergingAlgorithm::GetExtremalCoordinatesLowerX(const Cluster *const pCluster) 
{


  CartesianVector lowerXCoordinate(0,0,0);
  CartesianVector higherXCoordinate(0,0,0);

  GetExtremalCoordinates(pCluster, lowerXCoordinate, higherXCoordinate);

  return lowerXCoordinate;

}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector HitWidthClusterMergingAlgorithm::GetExtremalCoordinatesHigherX(const Cluster *const pCluster) 
{


  CartesianVector lowerXCoordinate(0,0,0);
  CartesianVector higherXCoordinate(0,0,0);

  GetExtremalCoordinates(pCluster, lowerXCoordinate, higherXCoordinate);

  return higherXCoordinate;

}


//------------------------------------------------------------------------------------------------------------------------------------------

//SHOULD REALLY BE IN A HELPER CLASS
void HitWidthClusterMergingAlgorithm::GetExtremalCoordinates(const Cluster *const pCluster, CartesianVector &lowerXCoordinate, CartesianVector &higherXCoordinate)
{
    return GetExtremalCoordinates(pCluster->GetOrderedCaloHitList(), lowerXCoordinate, higherXCoordinate);
}

//------------------------------------------------------------------------------------------------------------------------------------------


//SHOULD REALLY BE IN A HELPER CLASS
void HitWidthClusterMergingAlgorithm::GetExtremalCoordinates(const OrderedCaloHitList &orderedCaloHitList, CartesianVector &lowerXCoordinate, CartesianVector &higherXCoordinate)
{
    if (orderedCaloHitList.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    CartesianPointVector coordinateVector;

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        for(CaloHitList::const_iterator hitIter = iter->second->begin(); hitIter != iter->second->end(); ++hitIter) 
        {
            const CaloHit *const pCaloHit = *hitIter;

	    const float hitWidth = pCaloHit->GetCellSize1();
            const unsigned int numberOfConstituentHits = floor(hitWidth/pCaloHit->GetCellSize0()) + 1;
            const float constituentHitWidth = hitWidth/static_cast<float>(numberOfConstituentHits);
	 
	    // start at end of cluster with the lowest x value
            float xPositionAlongHit(pCaloHit->GetPositionVector().GetX() - (hitWidth/2));
	    for(unsigned int i(0); i < numberOfConstituentHits; ++i) 
	    {
                i == 0 ? xPositionAlongHit += constituentHitWidth/2 : xPositionAlongHit += constituentHitWidth;
	        CartesianVector consituentHitPosition(xPositionAlongHit, 0, pCaloHit->GetPositionVector().GetZ());
                coordinateVector.push_back(consituentHitPosition);
	    }
	}

    }

    std::sort(coordinateVector.begin(), coordinateVector.end(), LArClusterHelper::SortCoordinatesByPosition);
    return GetExtremalCoordinates(coordinateVector, lowerXCoordinate, higherXCoordinate);
}

//------------------------------------------------------------------------------------------------------------------------------------------



//SHOULD REALLY BE IN A HELPER CLASS
void HitWidthClusterMergingAlgorithm::GetExtremalCoordinates(const CartesianPointVector &coordinateVector, CartesianVector &lowerXCoordinate, CartesianVector &higherXCoordinate)
{
    if (coordinateVector.empty())
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

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


//SHOULD REALLY BE IN A HELPER CLASS
void HitWidthClusterMergingAlgorithm::GetExtremalCoordinatesX(const Cluster *const pCluster, float &minX, float &maxX)
{

    minX = +std::numeric_limits<float>::max();
    maxX = -std::numeric_limits<float>::max();

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(); hitIter  != iter->second->end(); ++hitIter)
        {

	    const CaloHit *const hit = (*hitIter);
	    float hitWidth = hit->GetCellSize1();

            unsigned int numberOfConstituentHits = floor(hitWidth/hit->GetCellSize0()) + 1;
            float constituentHitWidth = hitWidth/numberOfConstituentHits;

	    // lower has lower x value
	    float lowerHitEdge(hit->GetPositionVector().GetX() - (hitWidth/2.f));
	    float upperHitEdge(hit->GetPositionVector().GetX() + (hitWidth/2.f));

	    // First new point in the cluster
            minX = std::min(lowerHitEdge + constituentHitWidth/2.f, minX);
	    maxX = std::max(upperHitEdge - constituentHitWidth/2.f, maxX);

        }
    }

    if (maxX < minX)
        throw pandora::StatusCodeException(STATUS_CODE_FAILURE);

}

//------------------------------------------------------------------------------------------------------------------------------------------

//SHOULD REALLY BE IN A HELPER CLASS
bool HitWidthClusterMergingAlgorithm::SortByX(const Cluster *const pLhs, const Cluster *const pRhs)
{
    float lhsMinX(GetMinX(pLhs));
    float rhsMinX(GetMinX(pRhs));

    return (lhsMinX < rhsMinX);
}

//------------------------------------------------------------------------------------------------------------------------------------------

//SHOULD REALLY BE IN A HELPER CLASS
float HitWidthClusterMergingAlgorithm::GetMinX(const Cluster *const pCluster)
{

    float minX(+std::numeric_limits<float>::max());

    const OrderedCaloHitList *hits(&pCluster->GetOrderedCaloHitList());

    // find min x coordinate for each cluster
    for(OrderedCaloHitList::const_iterator iter = hits->begin(); iter !=  hits->end(); ++iter) 
    {
        for(CaloHitList::const_iterator hitIter = iter->second->begin(); hitIter != iter->second->end(); ++hitIter) 
        {
	    const CaloHit *const hit = (*hitIter);
	    float hitWidth = hit->GetCellSize1();

            unsigned int numberOfConstituentHits = floor(hitWidth/hit->GetCellSize0()) + 1;
            float constituentHitWidth = hitWidth/numberOfConstituentHits;

	    float lowHitEdge(hit->GetPositionVector().GetX() - (hitWidth/2));

	    // First new point in the cluster
            minX = std::min(lowHitEdge + constituentHitWidth/2, minX);
	}
    }
    return minX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HitWidthClusterMergingAlgorithm::SortByMaxX(const Cluster *const pLhs, const Cluster *const pRhs)
{

    CartesianVector lhsLowerXEdge(0,0,0);
    CartesianVector lhsUpperXEdge(0,0,0);

    CartesianVector rhsLowerXEdge(0,0,0);
    CartesianVector rhsUpperXEdge(0,0,0);

    HitWidthClusterMergingAlgorithm::GetExtremalCoordinates(pLhs, lhsLowerXEdge, lhsUpperXEdge);
    HitWidthClusterMergingAlgorithm::GetExtremalCoordinates(pRhs, rhsLowerXEdge, rhsUpperXEdge);
  
    return (lhsUpperXEdge.GetX() < rhsUpperXEdge.GetX());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float HitWidthClusterMergingAlgorithm::GetTotalClusterWeight(const Cluster *const pCluster)
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

void HitWidthClusterMergingAlgorithm::FillPositionToWeightMap(const pandora::Cluster *const pCluster, ClusterFitParameters::ClusterPositionToWeightMap &clusterPositionToWeightMap) const
{

    OrderedCaloHitList orderedCaloHitList = pCluster->GetOrderedCaloHitList();

    for(OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter !=  orderedCaloHitList.end(); ++iter) 
    {
        for(CaloHitList::const_iterator hitIter = iter->second->begin(); hitIter != iter->second->end(); ++hitIter) 
        {
	    const CaloHit *const hit = (*hitIter);

	    float hitWidth = hit->GetCellSize1();
	    float hitWeight = hitWidth;

            unsigned int numberOfConstituentHits = floor(hitWidth/hit->GetCellSize0()) + 1;
            float constituentHitWidth = hitWidth/numberOfConstituentHits;
	    float constituentHitWeight = hitWeight/numberOfConstituentHits;

	    float xPositionAlongHit(hit->GetPositionVector().GetX() - (hitWidth/2));
	    for(unsigned int i(0); i < numberOfConstituentHits; ++i) 
	    {
                i == 0 ? xPositionAlongHit += constituentHitWidth/2 : xPositionAlongHit += constituentHitWidth;
	         clusterPositionToWeightMap.insert(std::pair(CartesianVector(xPositionAlongHit, 0, hit->GetPositionVector().GetZ()), constituentHitWeight));
	    }
	}
    }

    return;

}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitWidthClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
 
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterWeight", m_minClusterWeight));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxXMergeDistance", m_maxXMergeDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxZMergeDistance", m_maxZMergeDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxMergeCosOpeningAngle", m_maxMergeCosOpeningAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FittingSampleWeight", m_fittingSampleWeight));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}
  
} // namespace lar_content



/*
  StatusCode HitWidthClusterMergingAlgorithm::Run(){

    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));

    ClusterVector clusterVector;
    GetListOfCleanClusters(pClusterList, clusterVector);

    
    // Fill cluster fit map 
    ClusterToFitParametersMap clusterToFitParametersMap;

    for(const Cluster *const pCluster : clusterVector) 
    {
      clusterToFitParametersMap.insert(std::pair(pCluster, ClusterFitParameters(pCluster)));
    }


    PandoraMonitoringApi::Create(this->GetPandora());
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);

    for(auto entry: clusterToFitParametersMap){

      ClusterList aCluster;
      aCluster.push_back(entry.first);

      PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &aCluster, "A Cluster", BLACK);

      std::cout << "Current: " << std::endl;

      for(auto iter = entry.second.m_currentClusterPositionToWeightMap.begin(); iter != entry.second.m_currentClusterPositionToWeightMap.end(); ++iter) {

	const CartesianVector pos = iter->first;

	PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &pos, "POINT", RED, 2);
	PandoraMonitoringApi::Pause(this->GetPandora());

      }

      PandoraMonitoringApi::Pause(this->GetPandora());


      std::cout << "Test: " << std::endl;

      for(auto iter = entry.second.m_testClusterPositionToWeightMap.begin(); iter != entry.second.m_testClusterPositionToWeightMap.end(); ++iter) {

	const CartesianVector pos = iter->first;

	PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &pos, "POINT", BLUE, 2);
	PandoraMonitoringApi::Pause(this->GetPandora());

      }

      PandoraMonitoringApi::Pause(this->GetPandora());

    }


   
    return STATUS_CODE_SUCCESS;

  } 
*/
//------------------------------------------------------------------------------------------------------------------------------------------

/*
StatusCode HitWidthClusterMergingAlgorithm::Run()
{

    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));

    ClusterVector clusterVector;
    GetListOfCleanClusters(pClusterList, clusterVector);
    
    PandoraMonitoringApi::Create(this->GetPandora());
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);

    //ClusterList consideredClusters(clusterVector.begin(), clusterVector.end());
    //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &consideredClusters, "Considered Clusters", BLACK);

    // Fill cluster fit map 
    ClusterToFitParametersMap clusterToFitParametersMap;
    for(const Cluster *const pCluster : clusterVector) 
    {
        if(!clusterToFitParametersMap.insert(std::pair(pCluster, ClusterFitParameters(pCluster))).second) {
	    std::cout << "Cluster fit not added to map!" << std::endl;
	    throw;
	}
    }



    for(ClusterVector::const_iterator iterCurrentCluster = clusterVector.begin(); iterCurrentCluster != clusterVector.end(); ++iterCurrentCluster)
    {
      
      const Cluster *const pCurrentCluster = *iterCurrentCluster;
      const ClusterFitParameters currentFitParameters(clusterToFitParametersMap.at(pCurrentCluster));
	
     ClusterList currentCluster;
     currentCluster.push_back(pCurrentCluster);
     PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentCluster, "Current Cluster", BLACK);

     const CartesianVector currentClusterDirection = GetClusterEndDirection(currentFitParameters.m_currentClusterPositionToWeightMap, currentFitParameters.m_numCaloHits);

      for (ClusterVector::const_iterator iterTestCluster = iterCurrentCluster; iterTestCluster != clusterVector.end(); ++iterTestCluster)
      {
	    
           if (iterCurrentCluster == iterTestCluster)
               continue;


	   const Cluster *const pTestCluster = *iterTestCluster;
           const ClusterFitParameters testFitParameters = clusterToFitParametersMap.at(pTestCluster);

           ClusterList testCluster;
           testCluster.push_back(pTestCluster);
           PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testCluster, "Test Cluster", BLUE);
	         
           const CartesianVector testClusterDirection = GetClusterEndDirection(testFitParameters.m_testClusterPositionToWeightMap, testFitParameters.m_numCaloHits);


	   std::cout << "Opening Angle: " << std::fabs(currentClusterDirection.GetCosOpeningAngle(testClusterDirection)) << std::endl;

	   PandoraMonitoringApi::Pause(this->GetPandora());

        }

      std::cout << "ALL CLUSTERS VIEWED" << std::endl;
      PandoraMonitoringApi::ViewEvent(this->GetPandora());

    }


    return STATUS_CODE_SUCCESS;
    
}
*/


//------------------------------------------------------------------------------------------------------------------------------------------
/*
StatusCode HitWidthClusterMergingAlgorithm::Run()
{

    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));

    ClusterVector clusterVector;
    GetListOfCleanClusters(pClusterList, clusterVector);
    
    PandoraMonitoringApi::Create(this->GetPandora());
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);

    //ClusterList consideredClusters(clusterVector.begin(), clusterVector.end());
    //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &consideredClusters, "Considered Clusters", BLACK);

    ClusterAssociationMap aMap;
    this->TestPopulateClusterAssociationMap(clusterVector, aMap);

    return STATUS_CODE_SUCCESS;
    
}
*/
