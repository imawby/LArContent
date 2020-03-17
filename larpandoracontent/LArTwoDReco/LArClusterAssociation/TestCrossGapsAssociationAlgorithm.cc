/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterAssociation/TestCrossGapsAssociationAlgorithm.cc
 *
 *  @brief  Implementation of the test cross gaps association algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterAssociation/TestCrossGapsAssociationAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TestCrossGapsAssociationAlgorithm::TestCrossGapsAssociationAlgorithm() :
    m_minClusterHits(10),
    m_minClusterLayers(6),
    m_slidingFitWindow(20),
    m_maxSamplingPoints(1000),
    m_sampleStepSize(0.5f),
    m_maxUnmatchedSampleRun(8),
    m_maxOnClusterDistance(1.5f),
    m_minMatchedSamplingPoints(10),
    m_minMatchedSamplingFraction(0.5f),
    m_gapTolerance(0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TestCrossGapsAssociationAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    // ATTN May want to opt-out completely if no gap information available
    // if (PandoraContentApi::GetGeometry(*this)->GetDetectorGapList().empty())
    //     return;


    PandoraMonitoringApi::Create(this->GetPandora());
    PandoraMonitoringApi::SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f);

    for (const Cluster *const pCluster : *pClusterList)
    {
      if (LArHitWidthHelper::GetTotalClusterWeight(pCluster) < 5) {
	  ClusterList aCluster;
	  aCluster.push_back(pCluster);
	  //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &aCluster, "CaloHits", RED);
	  //PandoraMonitoringApi::Pause(this->GetPandora());
            continue;
        }

	/*
        if (1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() < m_minClusterLayers) {
	  ClusterList aCluster;
	  aCluster.push_back(pCluster);
	  PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &aCluster, "Layer", RED);
	  PandoraMonitoringApi::Pause(this->GetPandora());
            continue;
	}
	*/

        clusterVector.push_back(pCluster);
    }

    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByInnerLayer);


    //ClusterList cleanClusters(clusterVector.begin(), clusterVector.end());
    //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &cleanClusters, "Clean Clusters", BLACK);

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());

}

//------------------------------------------------------------------------------------------------------------------------------------------

void TestCrossGapsAssociationAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{
    TwoDSlidingFitResultMap slidingFitResultMap;
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));

    
    for (const Cluster *const pCluster : clusterVector)
    {
        try {(void) slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch)));}
        catch (StatusCodeException &) {}
    }
    
    /*
    for (const Cluster *const pCluster : clusterVector)
    {
	
        CartesianPointVector uniformConstituentHits(LArHitWidthHelper::GetUniformConstituentHits(pCluster,0.5));

        try 
        {
	    TwoDSlidingFitResult fitResult(&uniformConstituentHits, m_slidingFitWindow, slidingFitPitch);
	    fitResult.SetCluster(pCluster);
            (void) slidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, fitResult));
            
        }
        catch (StatusCodeException &) {}
    }
    */
    // ATTN This method assumes that clusters have been sorted by layer
    for (ClusterVector::const_iterator iterI = clusterVector.begin(), iterIEnd = clusterVector.end(); iterI != iterIEnd; ++iterI)
    {
        const Cluster *const pInnerCluster = *iterI;
        TwoDSlidingFitResultMap::const_iterator fitIterI = slidingFitResultMap.find(pInnerCluster);;

        if (slidingFitResultMap.end() == fitIterI)
            continue;

        for (ClusterVector::const_iterator iterJ = iterI, iterJEnd = clusterVector.end(); iterJ != iterJEnd; ++iterJ)
        {
            const Cluster *const pOuterCluster = *iterJ;

            if (pInnerCluster == pOuterCluster)
                continue;

            TwoDSlidingFitResultMap::const_iterator fitIterJ = slidingFitResultMap.find(pOuterCluster);

            if (slidingFitResultMap.end() == fitIterJ)
                continue;

            ClusterList innerClusterList;
	    innerClusterList.push_back(pInnerCluster);

            ClusterList outerClusterList;
	    outerClusterList.push_back(pOuterCluster);

            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &innerClusterList, "Inner Cluster", DARKGREEN);
            //PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &outerClusterList, "Outer Cluster", DARKRED);

            if (!this->AreClustersAssociated(fitIterI->second, fitIterJ->second))
                continue;

            clusterAssociationMap[pInnerCluster].m_forwardAssociations.insert(pOuterCluster);
            clusterAssociationMap[pOuterCluster].m_backwardAssociations.insert(pInnerCluster);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TestCrossGapsAssociationAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster,  const Cluster *const pTestCluster) const
{
    const unsigned int currentLayer(isForward ? pCurrentCluster->GetOuterPseudoLayer() : pCurrentCluster->GetInnerPseudoLayer());
    const unsigned int testLayer(isForward ? pTestCluster->GetOuterPseudoLayer() : pTestCluster->GetInnerPseudoLayer());

    if (isForward && ((testLayer > currentLayer) || ((testLayer == currentLayer) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster))))
        return true;

    if (!isForward && ((testLayer < currentLayer) || ((testLayer == currentLayer) && LArClusterHelper::SortByNHits(pTestCluster, pCurrentCluster))))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TestCrossGapsAssociationAlgorithm::AreClustersAssociated(const TwoDSlidingFitResult &innerFitResult, const TwoDSlidingFitResult &outerFitResult) const
{
    if (outerFitResult.GetCluster()->GetInnerPseudoLayer() < innerFitResult.GetCluster()->GetInnerPseudoLayer())
        throw pandora::StatusCodeException(STATUS_CODE_NOT_ALLOWED);
    /*
    //not great for transverse tracks
    if (outerFitResult.GetCluster()->GetInnerPseudoLayer() < innerFitResult.GetCluster()->GetOuterPseudoLayer())
        return false;
    */


    // find the two 'extremal' cluster points that are closest together
    CartesianVector innerMergePosition(0,0,0);
    CartesianVector outerMergePosition(0,0,0);
    CartesianVector innerMergeDirection(0,0,0);
    CartesianVector outerMergeDirection(0,0,0);

    float mergePointsSeparation(std::numeric_limits<float>::max());

    CartesianVector innerMinLayerPosition(innerFitResult.GetGlobalMinLayerPosition());
    CartesianVector innerMaxLayerPosition(innerFitResult.GetGlobalMaxLayerPosition());
    CartesianVector outerMinLayerPosition(outerFitResult.GetGlobalMinLayerPosition());
    CartesianVector outerMaxLayerPosition(outerFitResult.GetGlobalMaxLayerPosition());

    //CartesianVector innerMinLayerPosition(std::numeric_limits<float>::max(),0,0);
    //CartesianVector innerMaxLayerPosition(-std::numeric_limits<float>::max(),0,0);
    //CartesianVector outerMinLayerPosition(std::numeric_limits<float>::max(),0,0);
    //CartesianVector outerMaxLayerPosition(-std::numeric_limits<float>::max(),0,0);

    //CartesianPointVector innerClusterHits(LArHitWidthHelper::GetUniformConstituentHits(innerFitResult.GetCluster(),0.5));
    //CartesianPointVector outerClusterHits(LArHitWidthHelper::GetUniformConstituentHits(outerFitResult.GetCluster(),0.5));

    //float innerMaxL(-std::numeric_limits<float>::max());
    //float innerMinL(std::numeric_limits<float>::max());

    //float outerMaxL(-std::numeric_limits<float>::max());
    //float outerMinL(std::numeric_limits<float>::max());
    /*
    float positionL(0);
    float positionR(0);
    for(const CartesianVector &position : innerClusterHits)
    {
      innerFitResult.GetLocalPosition(position, positionL, positionR);

      innerMinL = std::min(positionL, innerMinL);
      innerMaxL = std::max(positionL, innerMaxL);
    }

    for(const CartesianVector &position : outerClusterHits)
    {
      outerFitResult.GetLocalPosition(position, positionL, positionR);

      outerMinL = std::min(positionL, outerMinL);
      outerMaxL = std::max(positionL, outerMaxL);
    }
    */



    //LArClusterHelper::GetExtremalCoordinates(LArHitWidthHelper::GetUniformConstituentHits(innerFitResult.GetCluster(),0.5),  innerMinLayerPosition, innerMaxLayerPosition);
    //LArClusterHelper::GetExtremalCoordinates(LArHitWidthHelper::GetUniformConstituentHits(outerFitResult.GetCluster(),0.5),  outerMinLayerPosition, outerMaxLayerPosition);
    

    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,  innerFitResult.GetGlobalFitPosition(innerMinL, innerMinLayerPosition));
    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,  innerFitResult.GetGlobalFitPosition(innerMaxL, innerMaxLayerPosition));
    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,  outerFitResult.GetGlobalFitPosition(outerMinL, outerMinLayerPosition));
    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,  outerFitResult.GetGlobalFitPosition(outerMaxL, outerMaxLayerPosition));
    

    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &innerMinLayerPosition, "innerMinLayerPosition", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &innerMaxLayerPosition, "innerMaxLayerPosition", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &outerMinLayerPosition, "outerMinLayerPosition", RED, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &outerMaxLayerPosition, "outerMaxLayerPosition", RED, 2);
    //PandoraMonitoringApi::ViewEvent(this->GetPandora());

    
    CartesianVector innerMinLayerDirection(innerFitResult.GetGlobalMinLayerDirection());
    CartesianVector innerMaxLayerDirection(innerFitResult.GetGlobalMaxLayerDirection());
    CartesianVector outerMinLayerDirection(outerFitResult.GetGlobalMinLayerDirection());
    CartesianVector outerMaxLayerDirection(outerFitResult.GetGlobalMaxLayerDirection());

    //CartesianVector innerMinLayerDirection(0,0,0);
    //CartesianVector innerMaxLayerDirection(0,0,0);
    //CartesianVector outerMinLayerDirection(0,0,0);
    //CartesianVector outerMaxLayerDirection(0,0,0);
    
    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,  innerFitResult.GetGlobalFitDirection(innerMinL, innerMinLayerDirection));
    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,  innerFitResult.GetGlobalFitDirection(innerMaxL, innerMaxLayerDirection));
    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,  outerFitResult.GetGlobalFitDirection(outerMinL, outerMinLayerDirection));
    //PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=,  outerFitResult.GetGlobalFitDirection(outerMaxL, outerMaxLayerDirection));
    
    if((innerMinLayerPosition - outerMinLayerPosition).GetMagnitude() < mergePointsSeparation)
    {
        innerMergePosition = innerMinLayerPosition;
        outerMergePosition = outerMinLayerPosition;
        innerMergeDirection = innerMinLayerDirection;
	outerMergeDirection = outerMinLayerDirection;
	mergePointsSeparation = (innerMinLayerPosition - outerMinLayerPosition).GetMagnitude();
    }

    if((innerMinLayerPosition - outerMaxLayerPosition).GetMagnitude() < mergePointsSeparation)
    {
        innerMergePosition = innerMinLayerPosition;
        outerMergePosition = outerMaxLayerPosition;
        innerMergeDirection = innerMinLayerDirection;
	outerMergeDirection = outerMaxLayerDirection;
	mergePointsSeparation = (innerMinLayerPosition - outerMaxLayerPosition).GetMagnitude();
    }

    if((innerMaxLayerPosition - outerMinLayerPosition).GetMagnitude() < mergePointsSeparation)
    {
        innerMergePosition = innerMaxLayerPosition;
        outerMergePosition = outerMinLayerPosition;
        innerMergeDirection = innerMaxLayerDirection;
	outerMergeDirection = outerMinLayerDirection*(-1);
	mergePointsSeparation = (innerMaxLayerPosition - outerMinLayerPosition).GetMagnitude();
    }

    if((innerMaxLayerPosition - outerMaxLayerPosition).GetMagnitude() < mergePointsSeparation)
    {
        innerMergePosition = innerMaxLayerPosition;
        outerMergePosition = outerMaxLayerPosition;
        innerMergeDirection = innerMaxLayerDirection;
	outerMergeDirection = outerMaxLayerDirection;
	mergePointsSeparation = (innerMaxLayerPosition - outerMaxLayerPosition).GetMagnitude();
    }

    //make sure directions are pointing at each other
    int innerScale(1);
    int outerScale(1);
    if((innerMergePosition - (outerMergeDirection*(1/mergePointsSeparation))).GetMagnitude() > mergePointsSeparation)
      outerMergeDirection*(-1);

    if((outerMergePosition - (innerMergeDirection*(1/mergePointsSeparation))).GetMagnitude() > mergePointsSeparation)
      innerMergeDirection*(-1);

    return (this->IsAssociated(innerMergePosition, innerMergeDirection * innerScale, outerFitResult) &&
        this->IsAssociated(outerMergePosition, outerMergeDirection * outerScale, innerFitResult));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TestCrossGapsAssociationAlgorithm::IsAssociated(const CartesianVector &startPosition, const CartesianVector &startDirection,
    const TwoDSlidingFitResult &targetFitResult) const
{

    const HitType hitType(LArClusterHelper::GetClusterHitType(targetFitResult.GetCluster()));

    unsigned int nSamplingPoints(0), nGapSamplingPoints(0), nMatchedSamplingPoints(0), nUnmatchedSampleRun(0);


    typedef std::vector<const LineGap*> LineGapVector;
    LineGapVector lineGapVector;


    for (const DetectorGap *const pDetectorGap : this->GetPandora().GetGeometry()->GetDetectorGapList())
    {
        const LineGap *pLineGap(nullptr);
	pLineGap = dynamic_cast<const LineGap*>(pDetectorGap);

	if(pLineGap)
	    lineGapVector.push_back(pLineGap);
    }


    const LineGap *const pLineGap(lineGapVector[0]); 
      

    float gapStartX(pLineGap->GetLineStartX());
    float gapEndX(pLineGap->GetLineEndX());
    //float gapWidth(std::fabs(gapStartX - gapEndX));

    float m_maxXDistanceFromGap(5);

    if(std::fabs(gapStartX - startPosition.GetX()) > m_maxXDistanceFromGap && std::fabs(gapEndX - startPosition.GetX()) > m_maxXDistanceFromGap) 
        return false;


    //keep the issue the other side?
    CartesianVector samplingStartPoint = startPosition;

    std::fabs(gapStartX - startPosition.GetX()) <  std::fabs(gapEndX - startPosition.GetX()) ? samplingStartPoint = CartesianVector(gapEndX, 0, startPosition.GetZ()) : samplingStartPoint = CartesianVector(gapStartX, 0, startPosition.GetZ());


    for (unsigned int iSample = 0; iSample < m_maxSamplingPoints; ++iSample)
    {
        ++nSamplingPoints;
        const CartesianVector samplingPoint(samplingStartPoint + startDirection * static_cast<float>(iSample) * m_sampleStepSize);

	PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &samplingPoint, "Sampling Point", BLACK, 2);
	//PandoraMonitoringApi::Pause(this->GetPandora());

        if (LArGeometryHelper::IsInGap(this->GetPandora(), samplingPoint, hitType, m_gapTolerance))
        {
	    std::cout << "I AM IN A GAP" << std::endl;
            ++nGapSamplingPoints;
            nUnmatchedSampleRun = 0; // ATTN Choose to also reset run when entering gap region
            continue;
        }

        if (this->IsNearCluster(samplingPoint, targetFitResult))
        {
            ++nMatchedSamplingPoints;
            nUnmatchedSampleRun = 0;
        }
        else if (++nUnmatchedSampleRun > m_maxUnmatchedSampleRun)
        {
	    std::cout << "REACHED UNMATCHED LIMIT" << std::endl;
            break;
        }
    }


    const float expectation((targetFitResult.GetGlobalMaxLayerPosition() - targetFitResult.GetGlobalMinLayerPosition()).GetMagnitude() / m_sampleStepSize);
    const float matchedSamplingFraction(expectation > 0.f ? static_cast<float>(nMatchedSamplingPoints) / expectation : 0.f);


    std::cout << "Matched Sampling Points: " << nMatchedSamplingPoints << std::endl;
    std::cout << "Matched Sampling Fraction: " << matchedSamplingFraction << std::endl;

    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    if ((nMatchedSamplingPoints > m_minMatchedSamplingPoints) || (matchedSamplingFraction > m_minMatchedSamplingFraction))
        return true;

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TestCrossGapsAssociationAlgorithm::IsNearCluster(const CartesianVector &samplingPoint, const TwoDSlidingFitResult &targetFitResult) const
{
    float rL(std::numeric_limits<float>::max()), rT(std::numeric_limits<float>::max());
    targetFitResult.GetLocalPosition(samplingPoint, rL, rT);

    CartesianVector fitPosition(0.f, 0.f, 0.f);

    if (STATUS_CODE_SUCCESS == targetFitResult.GetGlobalFitPosition(rL, fitPosition))
    {
        if ((fitPosition - samplingPoint).GetMagnitudeSquared() < m_maxOnClusterDistance * m_maxOnClusterDistance)
            return true;
    }

    CartesianVector fitPositionAtX(0.f, 0.f, 0.f);

    if (STATUS_CODE_SUCCESS == targetFitResult.GetGlobalFitPositionAtX(samplingPoint.GetX(), fitPositionAtX))
    {
        if ((fitPositionAtX - samplingPoint).GetMagnitudeSquared() < m_maxOnClusterDistance * m_maxOnClusterDistance)
            return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------


StatusCode TestCrossGapsAssociationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterHits", m_minClusterHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLayers", m_minClusterLayers));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxSamplingPoints", m_maxSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SampleStepSize", m_sampleStepSize));

    if (m_sampleStepSize < std::numeric_limits<float>::epsilon())
    {
        std::cout << "TestCrossGapsAssociationAlgorithm: Invalid value for SampleStepSize " << m_sampleStepSize << std::endl;
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxUnmatchedSampleRun", m_maxUnmatchedSampleRun));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxOnClusterDistance", m_maxOnClusterDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingPoints", m_minMatchedSamplingPoints));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMatchedSamplingFraction", m_minMatchedSamplingFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "GapTolerance", m_gapTolerance));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
