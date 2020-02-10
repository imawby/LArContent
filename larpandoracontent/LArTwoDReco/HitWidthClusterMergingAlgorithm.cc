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

HitWidthClusterMergingAlgorithm::HitWidthClusterMergingAlgorithm() :
  m_clusterListName(),
  m_minClusterWeight(0),  //ATTN - THIS HAS BEEN SET IN THE MAC XML FILE
  m_maxXMergeDistance(0), //ATTN - THIS HAS BEEN SET IN THE MAC XML FILE
  m_maxZMergeDistance(0),  //ATTN - THIS HAS BEEN SET IN THE MAC XML FILE
  m_maxMergeCosOpeningAngle(0) //ATTN - THIS HAS BEEN TUNED IN THE MAC XML FILE
{
}

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

    ClusterList consideredClusters(clusterVector.begin(), clusterVector.end());
    
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &consideredClusters, "Considered Clusters", BLACK);
    PandoraMonitoringApi::Pause(this->GetPandora());

    ClusterAssociationMap aMap;
    this->TestPopulateClusterAssociationMap(clusterVector, aMap);

    return STATUS_CODE_SUCCESS;
    
}
  
  */
void HitWidthClusterMergingAlgorithm::TestPopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{

    // ATTN This method assumes that clusters have been sorted by x position (low x -> high x)
    for (ClusterVector::const_iterator iterCurrentCluster = clusterVector.begin(); iterCurrentCluster != clusterVector.end(); ++iterCurrentCluster)
    {

        const Cluster *const pCurrentCluster = *iterCurrentCluster;

        ClusterList currentClusterList;
        currentClusterList.push_back(pCurrentCluster);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentClusterList, "CURRENT", BLUE);

        CartesianVector currentLowerXEdge(0,0,0);
        CartesianVector currentUpperXEdge(0,0,0);
            
        GetExtremalCoordinates(pCurrentCluster, currentLowerXEdge, currentUpperXEdge);
	
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
            const Cluster *const pTestCluster = *iterTestCluster;
	    
            if (iterCurrentCluster == iterTestCluster)
                continue;
	    
            if (!this->TestAreClustersAssociated(pCurrentCluster, pTestCluster))
	        continue;

            clusterAssociationMap[pCurrentCluster].m_forwardAssociations.insert(pTestCluster);
            clusterAssociationMap[pTestCluster].m_backwardAssociations.insert(pCurrentCluster);
        }

	PandoraMonitoringApi::ViewEvent(this->GetPandora()); 
	
    }
}



  bool HitWidthClusterMergingAlgorithm::TestAreClustersAssociated(const Cluster *const pCurrentCluster, const Cluster *const pTestCluster) const
{

    //assumed that moving to higher x (when tied moving to higher z)
    CartesianVector currentLowerXEdge(0,0,0);
    CartesianVector currentUpperXEdge(0,0,0);
    
    CartesianVector testLowerXEdge(0,0,0);
    CartesianVector testUpperXEdge(0,0,0);

    GetExtremalCoordinates(pCurrentCluster, currentLowerXEdge, currentUpperXEdge);
    GetExtremalCoordinates(pTestCluster, testLowerXEdge, testUpperXEdge);

    ClusterList testClusterList;
    testClusterList.push_back(pTestCluster);

    
    CartesianVector currentClusterDirection(0, 0, 0);
    CartesianVector testClusterDirection(0, 0, 0);

    GetClusterDirection(pCurrentCluster, currentClusterDirection);
    GetClusterDirection(pTestCluster, testClusterDirection);
    
    std::string stringTag = "TEST - " + std::to_string(fabs(currentClusterDirection.GetCosOpeningAngle(testClusterDirection)));

    if(testLowerXEdge.GetX() > (currentUpperXEdge.GetX() + m_maxXMergeDistance)) {  //|| testLowerXEdge.GetX() < (currentUpperXEdge.GetX() - m_maxXMergeDistance)) {
      PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, "TEST - Outside X", RED);
      PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testLowerXEdge, "TEST", RED, 2);
      return false;
    }

    if(testLowerXEdge.GetZ() > (currentUpperXEdge.GetZ() + m_maxZMergeDistance) || testLowerXEdge.GetZ() < (currentUpperXEdge.GetZ() - m_maxZMergeDistance)) {
      PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, "TEST - Outside Z", RED);
      PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testLowerXEdge, "TEST", RED, 2);
      return false;
    }

    if(fabs(currentClusterDirection.GetCosOpeningAngle(testClusterDirection)) < m_maxMergeCosOpeningAngle) {
      PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, stringTag, RED);
      PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testLowerXEdge, "TEST", RED, 2);
      return false;
    }
    
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, stringTag, GREEN);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testLowerXEdge, "TEST", GREEN, 2);
    return true;

}

  
/*
StatusCode HitWidthClusterMergingAlgorithm::Run()
{
  
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, pClusterList));

    ClusterVector clusterVector;
    GetListOfCleanClusters(pClusterList, clusterVector);

   
    PandoraMonitoringApi::Create(this->GetPandora());
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);

    ClusterList consideredClusters(clusterVector.begin(), clusterVector.end());
    
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &consideredClusters, "Considered Clusters", BLACK);
    PandoraMonitoringApi::Pause(this->GetPandora());


    
    ClusterAssociationMap aMap;
    this->PopulateClusterAssociationMap(clusterVector, aMap);
    
    for(const Cluster *const pCluster : clusterVector) 
    {
      
      if (pCluster->GetNCaloHits() < m_minClusterHits)
	 continue;
        
      //ClusterFitResult clusterFit;
      //ClusterFitHelper::FitFullCluster(pCluster, clusterFit);

      //CartesianVector clusterFitDirection = clusterFit.GetDirection();
      //CartesianVector clusterFitIntercept = clusterFit.GetIntercept();
      //float clusterFitGradient = (clusterFitDirection.GetZ()/clusterFitDirection.GetX());
      //CartesianVector clusterFitZIntercept(0, 0, clusterFitIntercept.GetZ() - clusterFitGradient*clusterFitIntercept.GetX());
      //CartesianVector clusterFitEndpoint(-400, 0, clusterFitGradient*(-400) + clusterFitZIntercept.GetZ());

      
      //std::cout << "Intercept: "  << clusterFitIntercept << std::endl;
      //std::cout << "Endpoint: " << clusterDirection << std::endl;
      
      
      //float gradient(0);
      //float intercept(0);
      //float chiSquared(0);

      //GetWeightedGradient(pCluster, gradient, intercept, chiSquared);
      //CartesianVector interceptPoint(0, 0, intercept);
      //CartesianVector endpoint(-400, 0, gradient*(-400) + intercept);
      

      //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &clusterFitZIntercept, &clusterFitEndpoint, "Cluster Helper Fit", RED, 2, 2);
      //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &interceptPoint, &endpoint, "Least Squares Fit", BLUE, 2, 2);
      //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &singleCluster, "A Cluster", BLACK);
      //std::cout << "ChiSquared: " << chiSquared << std::endl;
      //PandoraMonitoringApi::Pause(this->GetPandora());

      //CartesianVector LSClusterFitDirection(0,0,0);
      //float LSIntercept(0);
      //GetClusterDirection(pCluster, LSClusterFitDirection);
      //GetClusterZIntercept(pCluster, LSIntercept);

      //float clusterFitGradient = (LSClusterFitDirection.GetZ()/LSClusterFitDirection.GetX());
      //CartesianVector LSStartPoint(0, 0, LSIntercept);
      //CartesianVector LSEndPoint(-400, 0, clusterFitGradient*(-400) + LSIntercept);

      //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &LSStartPoint, &LSEndPoint, "CHECK", RED, 2, 2);

      //PandoraMonitoringApi::Pause(this->GetPandora());
      
    }
    
    
      
    
    for(ClusterVector::const_iterator currentIter = clusterVector.begin(); currentIter != clusterVector.end(); ++currentIter) {

      const Cluster *const currentCluster(*currentIter);

      ClusterList currentClusterList;
      currentClusterList.push_back(currentCluster);

      PandoraMonitoringApi::ViewEvent(this->GetPandora());
      PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentClusterList, "Current Cluster", BLACK);

      for(ClusterVector::const_iterator testIter = clusterVector.begin(); testIter != (clusterVector.end()); ++testIter) {

          if(currentIter == testIter)
	      continue;

	  const Cluster *const testCluster(*testIter);

          ClusterList testClusterList;
          testClusterList.push_back(testCluster);

          //bool areAssociated(this->AreClustersAssociated(currentCluster, testCluster));
	  bool isExtremalCluster(this->IsExtremalCluster(true, currentCluster, testCluster));
	  
	  //areAssociated? PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, "Test Cluster", DARKGREEN) : PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, "Test Cluster", DARKRED);

	  isExtremalCluster? PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, "Test Cluster", DARKGREEN) : PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, "Test Cluster", DARKRED);

      }

      PandoraMonitoringApi::Pause(this->GetPandora());
    }
   

    //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), pClusterList, "All Clusters", BLACK);
    //PandoraMonitoringApi::Pause(this->GetPandora());



    return STATUS_CODE_SUCCESS;
  
}
*/

void HitWidthClusterMergingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{

  //ignore clusters that aren't transverse??

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

    // ATTN This method assumes that clusters have been sorted by x position (low x -> high x)
    for (ClusterVector::const_iterator iterCurrentCluster = clusterVector.begin(); iterCurrentCluster != clusterVector.end(); ++iterCurrentCluster)
    {

        const Cluster *const pCurrentCluster = *iterCurrentCluster;

        //ClusterList currentClusterList;
        //currentClusterList.push_back(pCurrentCluster);
        //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &currentClusterList, "CURRENT", BLUE);

        CartesianVector currentLowerXEdge(0,0,0);
        CartesianVector currentUpperXEdge(0,0,0);
            
        GetExtremalCoordinates(pCurrentCluster, currentLowerXEdge, currentUpperXEdge);
	
        //CartesianVector mergeXDist(currentUpperXEdge.GetX() + m_maxXMergeDistance, 0, currentUpperXEdge.GetZ() - m_maxZMergeDistance);
        //CartesianVector mergeZDist1(currentUpperXEdge.GetX(), 0, currentUpperXEdge.GetZ() + m_maxZMergeDistance);
        //CartesianVector mergeZDist2(currentUpperXEdge.GetX(), 0, currentUpperXEdge.GetZ() - m_maxZMergeDistance);
        //CartesianVector farPoint(currentUpperXEdge.GetX() + m_maxXMergeDistance, 0, currentUpperXEdge.GetZ() + m_maxZMergeDistance);

        //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &mergeZDist2, &mergeZDist1, "1", RED, 2, 2);
        //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &mergeZDist2, &mergeXDist, "2", RED, 2, 2);
        //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &farPoint, &mergeZDist1, "3", RED, 2, 2);
        //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &farPoint, &mergeXDist, "4", RED, 2, 2);

	//PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &currentUpperXEdge, "Upper Coordinate", BLUE, 2);
	
        for (ClusterVector::const_iterator iterTestCluster = iterCurrentCluster; iterTestCluster != clusterVector.end(); ++iterTestCluster)
        {
            const Cluster *const pTestCluster = *iterTestCluster;
	    
            if (iterCurrentCluster == iterTestCluster)
                continue;
	    
            if (!this->AreClustersAssociated(pCurrentCluster, pTestCluster))
	        continue;

            clusterAssociationMap[pCurrentCluster].m_forwardAssociations.insert(pTestCluster);
            clusterAssociationMap[pTestCluster].m_backwardAssociations.insert(pCurrentCluster);
        }

	//PandoraMonitoringApi::ViewEvent(this->GetPandora()); 
	
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HitWidthClusterMergingAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster,  const Cluster *const pTestCluster) const
{
  
    float currentMinX(0), currentMaxX(0);
    GetExtremalCoordinatesX(pCurrentCluster, currentMinX, currentMaxX);

    float testMinX(0.f), testMaxX(0.f);
    GetExtremalCoordinatesX(pTestCluster, testMinX, testMaxX);

    if (isForward)
    {
        if (std::fabs(testMaxX - currentMaxX) > std::numeric_limits<float>::epsilon())
            return (testMaxX > currentMaxX);
    }
    else
    {
        if (std::fabs(testMinX - currentMaxX) > std::numeric_limits<float>::epsilon())
            return (testMinX < currentMinX);
    }
  
    return LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster);

}


//------------------------------------------------------------------------------------------------------------------------------------------


  bool HitWidthClusterMergingAlgorithm::AreClustersAssociated(const Cluster *const pCurrentCluster, const Cluster *const pTestCluster) const
{

    //assumed that moving to higher x (when tied moving to higher z)
    CartesianVector currentLowerXEdge(0,0,0);
    CartesianVector currentUpperXEdge(0,0,0);
    
    CartesianVector testLowerXEdge(0,0,0);
    CartesianVector testUpperXEdge(0,0,0);

    GetExtremalCoordinates(pCurrentCluster, currentLowerXEdge, currentUpperXEdge);
    GetExtremalCoordinates(pTestCluster, testLowerXEdge, testUpperXEdge);

    //ClusterList testClusterList;
    //testClusterList.push_back(pTestCluster);

    
    CartesianVector currentClusterDirection(0, 0, 0);
    CartesianVector testClusterDirection(0, 0, 0);

    GetClusterDirection(pCurrentCluster, currentClusterDirection);
    GetClusterDirection(pTestCluster, testClusterDirection);
    
    //std::string stringTag = "TEST - " + std::to_string(fabs(currentClusterDirection.GetCosOpeningAngle(testClusterDirection)));

    if(testLowerXEdge.GetX() > (currentUpperXEdge.GetX() + m_maxXMergeDistance)) {  //|| testLowerXEdge.GetX() < (currentUpperXEdge.GetX() - m_maxXMergeDistance)) {
      //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, "TEST - Outside X", RED);
      //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testLowerXEdge, "TEST", RED, 2);
      return false;
    }

    if(testLowerXEdge.GetZ() > (currentUpperXEdge.GetZ() + m_maxZMergeDistance) || testLowerXEdge.GetZ() < (currentUpperXEdge.GetZ() - m_maxZMergeDistance)) {
      //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, "TEST - Outside Z", RED);
      //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testLowerXEdge, "TEST", RED, 2);
      return false;
    }

    if(fabs(currentClusterDirection.GetCosOpeningAngle(testClusterDirection)) < m_maxMergeCosOpeningAngle) {
      //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, stringTag, RED);
      //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testLowerXEdge, "TEST", RED, 2);
      return false;
    }
    
    //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &testClusterList, stringTag, GREEN);
    //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &testLowerXEdge, "TEST", GREEN, 2);
    return true;

}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetClusterDirection(const Cluster *const pCluster, CartesianVector &direction) const 
{

    float gradient(0);
    float intercept(0);
    float chiSquared(0);
    
    GetWeightedGradient(pCluster, gradient, intercept, chiSquared);

    float directionX = pow(1 + pow(gradient, 2), -0.5);
    float directionZ = gradient * directionX;
    
    direction = CartesianVector(directionX, 0, directionZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetClusterZIntercept(const Cluster *const pCluster, float &intercept) const
{
    float gradient(0);
    float chiSquared(0);
    
    GetWeightedGradient(pCluster, gradient, intercept, chiSquared);

    return;

}

//------------------------------------------------------------------------------------------------------------------------------------------


void HitWidthClusterMergingAlgorithm::GetWeightedGradient(const Cluster *const pCluster, float &gradient, float &intercept, float &chiSquared) const
{

    OrderedCaloHitList orderedCaloHitList = pCluster->GetOrderedCaloHitList();

    float weightSum(0);
    float weightedXSum(0);
    float weightedZSum(0);

    // Loop through hits in cluster to find the weighted x sum and weighted z sum of hits
    for(OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter !=  orderedCaloHitList.end(); ++iter) 
    {
        for(CaloHitList::const_iterator hitIter = iter->second->begin(); hitIter != iter->second->end(); ++hitIter) 
        {

	    const CaloHit *const hit = (*hitIter);

	    float hitWidth = hit->GetCellSize1();

	    float hitWeight = hitWidth;

	    //unsigned int numberOfConstituentHits = 1;
            unsigned int numberOfConstituentHits = floor(hitWidth/hit->GetCellSize0()) + 1;
            float constituentHitWidth = hitWidth/numberOfConstituentHits;
	    float constituentHitWeight = hitWeight/numberOfConstituentHits;

	    float xPositionAlongHit(hit->GetPositionVector().GetX() - (hitWidth/2));
	    for(unsigned int i(0); i < numberOfConstituentHits; ++i) 
	    {
                i == 0 ? xPositionAlongHit += constituentHitWidth/2 : xPositionAlongHit += constituentHitWidth;
		
		weightedXSum += xPositionAlongHit*constituentHitWeight;

		//CartesianVector positionVector(xPositionAlongHit, 0, hit->GetPositionVector().GetZ());
                //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &positionVector, "Constituent Hit Centre", DARKGREEN, 2);
	    }

	    weightedZSum += hit->GetPositionVector().GetZ() * hitWeight;
	    weightSum += hitWeight;
	}
    }

    float weightedXMean(weightedXSum/weightSum);
    float weightedZMean(weightedZSum/weightSum); 


    float numerator(0);
    float denominator(0);

    for(OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter !=  orderedCaloHitList.end(); ++iter) 
    {
        for(CaloHitList::const_iterator hitIter = iter->second->begin(); hitIter != iter->second->end(); ++hitIter) 
        {

	    const CaloHit *const hit = (*hitIter);

	    float hitWidth = hit->GetCellSize1();

	    float hitWeight = hitWidth;

	    //unsigned int numberOfConstituentHits = 1;
            unsigned int numberOfConstituentHits = floor(hitWidth/hit->GetCellSize0()) + 1;
            float constituentHitWidth = hitWidth/numberOfConstituentHits;
	    float constituentHitWeight = hitWeight/numberOfConstituentHits;

	    float xPositionAlongHit(hit->GetPositionVector().GetX() - (hitWidth/2));
	    for(unsigned int i(0); i < numberOfConstituentHits; ++i) 
	    {
                i == 0 ? xPositionAlongHit += constituentHitWidth/2 : xPositionAlongHit += constituentHitWidth;
		
		numerator += constituentHitWeight*(xPositionAlongHit - weightedXMean)*(hit->GetPositionVector().GetZ() - weightedZMean);
		denominator += constituentHitWeight*pow(xPositionAlongHit - weightedXMean, 2); 

		//CartesianVector positionVector(xPositionAlongHit, 0, hit->GetPositionVector().GetZ());
                //PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &positionVector, "Constituent Hit Centre", DARKGREEN, 2);
	    }
	}

	gradient = numerator/denominator;
        intercept = weightedZMean - gradient*weightedXMean;

    }


    chiSquared = 0;

    return;

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

    GetExtremalCoordinates(pLhs, lhsLowerXEdge, lhsUpperXEdge);
    GetExtremalCoordinates(pRhs, rhsLowerXEdge, rhsUpperXEdge);
  
    return (lhsUpperXEdge.GetX() < rhsUpperXEdge.GetX());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float HitWidthClusterMergingAlgorithm::GetTotalClusterWeight(const Cluster *const pCluster)
{
  
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    float hitWeight(0);
    
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

    return STATUS_CODE_SUCCESS;
}
  
} // namespace lar_content




/*



   for(OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter !=  orderedCaloHitList.end(); ++iter) 
    {
        for(CaloHitList::const_iterator hitIter = iter->second->begin(); hitIter != iter->second->end(); ++hitIter) 
        {

	    float hitWidth = (*hitIter)->GetCellSize1();
	    float hitWeight = hitWidth;

	    std::cout << "Hit Width: " << hitWidth << std::endl;
	    std::cout << "CellSize0: " << (*hitIter)->GetCellSize0() << std::endl;

            widthSum += hitWidth;

            sumX += ((*hitIter)->GetPositionVector().GetX() - (hitWidth/2))*(hitWeight/3);
	    sumX += ((*hitIter)->GetPositionVector().GetX())*(hitWeight/3);
            sumX += ((*hitIter)->GetPositionVector().GetX() + (hitWidth/2))*(hitWeight/3);

            sumZ += ((*hitIter)->GetPositionVector().GetZ())*hitWeight;
	}
    }

    float meanX(sumX/widthSum);
    float meanZ(sumZ/widthSum);

    float numerator(0);
    float denominator(0);

    for(OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter !=  orderedCaloHitList.end(); ++iter) 
    {
        for(CaloHitList::const_iterator hitIter = iter->second->begin(); hitIter != iter->second->end(); ++hitIter) 
        {
	    float hitWidth = (*hitIter)->GetCellSize1();
            float hitWeight = hitWidth;

	    numerator += (((*hitIter)->GetPositionVector().GetX() - (hitWidth/2)) - meanX) * ((*hitIter)->GetPositionVector().GetZ() - meanZ) * (hitWeight/3);
	    numerator += (((*hitIter)->GetPositionVector().GetX() + (hitWidth/2)) - meanX) * ((*hitIter)->GetPositionVector().GetZ() - meanZ) * (hitWeight/3);
	    numerator += (((*hitIter)->GetPositionVector().GetX()) - meanX) * ((*hitIter)->GetPositionVector().GetZ() - meanZ) * (hitWeight/3);

	    denominator += (hitWeight/3)*pow((*hitIter)->GetPositionVector().GetX() - meanX, 2);
	    denominator += (hitWeight/3)*pow((*hitIter)->GetPositionVector().GetX() - (hitWidth/2) - meanX, 2);
	    denominator += (hitWeight/3)*pow((*hitIter)->GetPositionVector().GetX() + (hitWidth/2) - meanX, 2);
	}
    }

    gradient = numerator/denominator;

    intercept = meanZ - gradient*meanX;


    for(OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter !=  orderedCaloHitList.end(); ++iter) 
    {
        for(CaloHitList::const_iterator hitIter = iter->second->begin(); hitIter != iter->second->end(); ++hitIter) 
        {
	    float hitWidth = (*hitIter)->GetCellSize1();
            float hitWeight = hitWidth;

	    chiSquared += (hitWeight/3) * pow((*hitIter)->GetPositionVector().GetZ() - intercept - gradient*((*hitIter)->GetPositionVector().GetX() + hitWidth/2), 2);
	    chiSquared += (hitWeight/3) * pow((*hitIter)->GetPositionVector().GetZ() - intercept - gradient*((*hitIter)->GetPositionVector().GetX() - hitWidth/2), 2);
	    chiSquared += (hitWeight/3) * pow((*hitIter)->GetPositionVector().GetZ() - intercept - gradient*((*hitIter)->GetPositionVector().GetX()), 2);

	}
    }



 */
