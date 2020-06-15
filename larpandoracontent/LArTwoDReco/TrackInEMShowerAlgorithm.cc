/**
 *  @file   TrackInEMShowerAlgorithm.cc
 *
 *  @brief  Implementation of the track in em shower algorithm class
 *
 *  $Log: $
 */
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArTwoDReco/TrackInEMShowerAlgorithm.h"

using namespace pandora;

namespace lar_content
{

TrackInEMShowerAlgorithm::ClusterAssociation::ClusterAssociation(const Cluster *const pUpstreamCluster, const Cluster *const pDownstreamCluster,const CartesianVector &upstreamMergePoint,
        const CartesianVector &upstreamMergeDirection, const CartesianVector &downstreamMergePoint, const CartesianVector &downstreamMergeDirection) :
    m_pUpstreamCluster(pUpstreamCluster),
    m_pDownstreamCluster(pDownstreamCluster),
    m_upstreamMergePoint(upstreamMergePoint),
    m_upstreamMergeDirection(upstreamMergeDirection),
    m_downstreamMergePoint(downstreamMergePoint),
    m_downstreamMergeDirection(downstreamMergeDirection),
    m_connectingLineDirection(0.f, 0.f, 0.f)
{
    const CartesianVector connectingLineDirection(m_downstreamMergePoint.GetX() - m_upstreamMergePoint.GetX(), 0.f, m_downstreamMergePoint.GetZ() - m_upstreamMergePoint.GetZ());
    m_connectingLineDirection = connectingLineDirection.GetUnitVector();
}

//------------------------------------------------------------------------------------------------------------------------------------------

TrackInEMShowerAlgorithm::ClusterAssociation::ClusterAssociation() :
    m_pUpstreamCluster(nullptr),
    m_pDownstreamCluster(nullptr),
    m_upstreamMergePoint(CartesianVector(0.f, 0.f, 0.f)),
    m_upstreamMergeDirection(CartesianVector(0.f, 0.f, 0.f)),
    m_downstreamMergePoint(CartesianVector(0.f, 0.f, 0.f)),
    m_downstreamMergeDirection(CartesianVector(0.f, 0.f, 0.f)),
    m_connectingLineDirection(CartesianVector(0.f, 0.f, 0.f))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------
    
TrackInEMShowerAlgorithm::TrackInEMShowerAlgorithm() :
    m_minCaloHits(25), 
    m_minSeparationDistance(27.f),
    m_maxPredictedMergePointOffset(5.f),
    m_slidingFitWindow(20),
    m_mergePointMinCosAngleDeviation(0.9995),
    m_minClusterLengthSum(75.f),
    m_minDirectionDeviationCosAngle(0.99),
    m_distanceFromLine(0.35),
    m_maxTrackGaps(2),
    m_lineSegmentLength(7.f),
    m_maxHitDistanceFromCluster(3.f),
    m_maxHitsToRecluster(4),
    m_minModificationFractionToRecluster(0.8)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
bool TrackInEMShowerAlgorithm::SortByDistanceAlongLine::operator() (const pandora::CaloHit *const pLhs, const pandora::CaloHit *const pRhs)
{
    const CartesianVector lhsDistanceVector(pLhs->GetPositionVector() - m_startPoint);
    const CartesianVector rhsDistanceVector(pRhs->GetPositionVector() - m_startPoint);

    const float lhsProjectedDistance(lhsDistanceVector.GetDotProduct(m_lineDirection));
    const float rhsProjectedDistance(rhsDistanceVector.GetDotProduct(m_lineDirection));

    return (lhsProjectedDistance < rhsProjectedDistance);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackInEMShowerAlgorithm::Run()
{
    const ClusterList *pClusterList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    const CaloHitList *pCaloHitList = NULL;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));
    
    ClusterVector clusterVector;
    this->SelectCleanClusters(pClusterList, clusterVector);

    TwoDSlidingFitResultMap microSlidingFitResultMap;
    TwoDSlidingFitResultMap macroSlidingFitResultMap;
    this->InitialiseSlidingFitResultMaps(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

    // Find longest length association and conitnue merging until no further merges can be made
    bool mergeMade(false);
    do
    {
        ClusterAssociation clusterAssociation;
        if(!this->FindBestClusterAssociation(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap, clusterAssociation))
            break;
        
        ClusterToCaloHitListMap clusterToCaloHitListMap;
        CaloHitVector extrapolatedCaloHitVector;
        this->GetExtrapolatedCaloHits(clusterAssociation, pClusterList, extrapolatedCaloHitVector, clusterToCaloHitListMap);
        
        if (extrapolatedCaloHitVector.empty())
            break;
        
        if (!this->IsTrackContinuous(clusterAssociation, extrapolatedCaloHitVector))
            break;

        ClusterList innerCluster({clusterAssociation.GetUpstreamCluster()}), outerCluster({clusterAssociation.GetDownstreamCluster()});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &innerCluster, "CLUSTER", BLACK);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &outerCluster, "CLUSTER", BLACK);
        
        this->RefineTracks(clusterAssociation, microSlidingFitResultMap, macroSlidingFitResultMap, clusterVector, extrapolatedCaloHitVector);

        for (const CaloHit *const pCaloHit : extrapolatedCaloHitVector)
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &hitPosition, "EXTRAPOLATED HIT", GREEN, 2);
        }

        
        this->MergeHits(clusterAssociation, clusterToCaloHitListMap, microSlidingFitResultMap, macroSlidingFitResultMap, clusterVector);

        this->UpdateSlidingFitResultMap(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

        mergeMade = true;
    } while (mergeMade);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::SelectCleanClusters(const ClusterList *pClusterList, ClusterVector &clusterVector) const
{
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster->GetNCaloHits() < m_minCaloHits)
            continue;

        clusterVector.push_back(pCluster);
    }
    
    std::sort(clusterVector.begin(), clusterVector.end(), LArClusterHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::InitialiseSlidingFitResultMaps(const ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
    TwoDSlidingFitResultMap &macroSlidingFitResultMap) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const unsigned int macroSlidingFitWindow(1000);
    
    for (const Cluster *const pCluster : clusterVector)
    {
        try
        {
            (void) microSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch)));
            (void) macroSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, macroSlidingFitWindow, slidingFitPitch)));
        }
        catch (StatusCodeException &) {}
    }
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::FindBestClusterAssociation(ClusterVector &clusterVector, const TwoDSlidingFitResultMap &microSlidingFitResultMap,
    const TwoDSlidingFitResultMap &macroSlidingFitResultMap, ClusterAssociation &clusterAssociation) const
{
    bool foundAssociation(false);
    float maxLength(0.f);

    // Find the associated cluster pair with the longest length sum
    for (ClusterVector::const_iterator currentIter = clusterVector.begin(); currentIter != clusterVector.end(); ++currentIter)
    {
        const Cluster *const pCurrentCluster(*currentIter);
        
        const TwoDSlidingFitResultMap::const_iterator currentMicroFitIter(microSlidingFitResultMap.find(pCurrentCluster));
        if (currentMicroFitIter == microSlidingFitResultMap.end())
            return false;

        const TwoDSlidingFitResultMap::const_iterator currentMacroFitIter(macroSlidingFitResultMap.find(pCurrentCluster));
        if (currentMacroFitIter == macroSlidingFitResultMap.end())
            return false;
   
        for (ClusterVector::const_iterator testIter = std::next(currentIter); testIter != clusterVector.end(); ++testIter)
        {
            const Cluster *const pTestCluster(*testIter);

            const float lengthSum(LArClusterHelper::GetLength(pCurrentCluster) + LArClusterHelper::GetLength(pTestCluster)); 
            if ((lengthSum < maxLength) || (lengthSum < m_minClusterLengthSum))
                continue;

            const TwoDSlidingFitResultMap::const_iterator testMicroFitIter(microSlidingFitResultMap.find(pTestCluster));
            if (testMicroFitIter == microSlidingFitResultMap.end())
                continue;

            const TwoDSlidingFitResultMap::const_iterator testMacroFitIter(macroSlidingFitResultMap.find(pTestCluster));
            if (testMacroFitIter == macroSlidingFitResultMap.end())
                continue;

            const bool isCurrentUpstream(LArClusterHelper::SortByPosition(pCurrentCluster, pTestCluster));

            // Ensure that clusters are not contained within one another
            const float currentMinLayerZ(currentMacroFitIter->second.GetGlobalMinLayerPosition().GetZ()), currentMaxLayerZ(currentMacroFitIter->second.GetGlobalMaxLayerPosition().GetZ());
            const float testMinLayerZ(testMacroFitIter->second.GetGlobalMinLayerPosition().GetZ()), testMaxLayerZ(testMacroFitIter->second.GetGlobalMaxLayerPosition().GetZ());
            
            if (((currentMinLayerZ > testMinLayerZ) && (currentMaxLayerZ < testMaxLayerZ)) || ((testMinLayerZ > currentMinLayerZ) && (testMaxLayerZ < currentMaxLayerZ)))
                continue;

            CartesianVector currentMergePoint(0.f, 0.f, 0.f), testMergePoint(0.f, 0.f, 0.f), currentMergeDirection(0.f, 0.f, 0.f), testMergeDirection(0.f, 0.f, 0.f);
            if (!this->GetClusterMergingCoordinates(currentMicroFitIter->second, currentMacroFitIter->second, testMacroFitIter->second, currentMergePoint, currentMergeDirection, isCurrentUpstream) ||
                !this->GetClusterMergingCoordinates(testMicroFitIter->second, testMacroFitIter->second, currentMacroFitIter->second, testMergePoint, testMergeDirection, !isCurrentUpstream))
            {
                continue;
            }

            if ((isCurrentUpstream && !this->AreClustersAssociated(currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection)) ||
                (!isCurrentUpstream && !this->AreClustersAssociated(testMergePoint, testMergeDirection, currentMergePoint, currentMergeDirection)))
            {
                continue;
            }

            foundAssociation = true;
            maxLength = lengthSum;

            if (isCurrentUpstream)
            {
                clusterAssociation = ClusterAssociation(pCurrentCluster, pTestCluster, currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection);
            }
            else
            {
                clusterAssociation = ClusterAssociation(pTestCluster, pCurrentCluster, testMergePoint, testMergeDirection, currentMergePoint, currentMergeDirection);
            }
        }
    }

    return foundAssociation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::GetClusterMergingCoordinates(const TwoDSlidingFitResult &currentMicroFitResult, const TwoDSlidingFitResult &currentMacroFitResult,
    const TwoDSlidingFitResult &associatedMacroFitResult, CartesianVector &currentMergePosition, CartesianVector &currentMergeDirection, const bool isUpstream) const
{
    CartesianVector currentAverageDirection(0.f, 0.f, 0.f);
    CartesianVector associatedAverageDirection(0.f, 0.f, 0.f);
    currentMacroFitResult.GetGlobalDirection(currentMacroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), currentAverageDirection);
    associatedMacroFitResult.GetGlobalDirection(associatedMacroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), associatedAverageDirection);    
    
    const LayerFitResultMap& currentMicroLayerFitResultMap(currentMicroFitResult.GetLayerFitResultMap());
    const unsigned int startLayer(isUpstream ? currentMicroFitResult.GetMaxLayer() : currentMicroFitResult.GetMinLayer());
    const unsigned int endLayer(isUpstream ? currentMicroFitResult.GetMinLayer() : currentMicroFitResult.GetMaxLayer());
    const unsigned int loopTerminationLayer(endLayer + (isUpstream ? -1 : 1));
    const int step(isUpstream ? -1 : 1);
    unsigned int goodPositionCount(0);
    unsigned int gradientStabilityHitWindow(std::ceil(currentMicroFitResult.GetCluster()->GetNCaloHits() * 0.1));
    
    for (unsigned int i = startLayer; i != loopTerminationLayer; i += step)
    {
        const auto microIter(currentMicroLayerFitResultMap.find(i));

        if (microIter == currentMicroLayerFitResultMap.end())
            continue;
        
        CartesianVector microDirection(0.f, 0.f, 0.f);
        currentMicroFitResult.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
    
        const float cosDirectionOpeningAngle(microDirection.GetCosOpeningAngle(associatedAverageDirection));
        if (cosDirectionOpeningAngle > m_mergePointMinCosAngleDeviation)
        {
            if (goodPositionCount == 0)
            {
                // so that direction vectors face one another
                // ISOBEL WHY ARE YOU USING THE AVERAGE DIRECTION???????
                currentMergeDirection = currentAverageDirection * (isUpstream ? 1.f : -1.f);
                currentMicroFitResult.GetGlobalFitPosition(microIter->second.GetL(), currentMergePosition);
            }
            
            ++goodPositionCount;
        }
        else
        {
            goodPositionCount = 0;
        }

        if (goodPositionCount > gradientStabilityHitWindow)
            break;
        
        // abort merging process if cannot find a stable region with gradient close to associated average 
        if (i == endLayer)
            return false;                                         
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::AreClustersAssociated(const CartesianVector &upstreamPoint, const CartesianVector &upstreamDirection, const CartesianVector &downstreamPoint,
    const CartesianVector &downstreamDirection) const
{
    if (downstreamPoint.GetZ() < upstreamPoint.GetZ())
        return false;

    // check that clusters are reasonably far away
    const float separationDistance(std::sqrt(upstreamPoint.GetDistanceSquared(downstreamPoint)));
    if (separationDistance < m_minSeparationDistance)
        return false;
    
    // check that opening angle is not too large
    // ATTN: cluster directions are pointing towards one another
    if (upstreamDirection.GetCosOpeningAngle(downstreamDirection * (-1.0)) < m_minDirectionDeviationCosAngle)
        return false;
    
    // check that fit allows you to get from upstream merge point to downstream point
    const CartesianVector predictedDownstreamPoint(upstreamPoint + (upstreamDirection * separationDistance));
    const float predictedDownstreamOffset((predictedDownstreamPoint - downstreamPoint).GetMagnitude());
    
    if (predictedDownstreamOffset > m_maxPredictedMergePointOffset)
        return false;

    // check that fit allows you to get from downstream point to upstream point
    const CartesianVector predictedUpstreamPoint(downstreamPoint + (downstreamDirection * separationDistance));
    const float predictedUpstreamOffset((predictedUpstreamPoint - upstreamPoint).GetMagnitude());

    if (predictedUpstreamOffset > m_maxPredictedMergePointOffset)
        return false;
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::GetExtrapolatedCaloHits(const ClusterAssociation &clusterAssociation, const ClusterList *const pClusterList,
    CaloHitVector &extrapolatedCaloHitVector, ClusterToCaloHitListMap &clusterToCaloHitListMap) const
{
    const CartesianVector &upstreamPoint(clusterAssociation.GetUpstreamMergePoint()), &downstreamPoint(clusterAssociation.GetDownstreamMergePoint());
    const float minX(std::min(upstreamPoint.GetX(), downstreamPoint.GetX())), maxX(std::max(upstreamPoint.GetX(), downstreamPoint.GetX()));

    // ATTN: consider hits from merging clusters as merging point may not be cluster end
    for (const Cluster *const pCluster : *pClusterList) 
    {
        OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList());
        for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
        {
            for (const CaloHit *const pCaloHit : *mapEntry.second) 
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

                if ((hitPosition.GetX() < minX) || (hitPosition.GetX() > maxX) || (hitPosition.GetZ() < upstreamPoint.GetZ()) || (hitPosition.GetZ() > downstreamPoint.GetZ()))
                    continue;
                
                const float distanceFromLine(clusterAssociation.GetConnectingLineDirection().GetCrossProduct(hitPosition - upstreamPoint).GetMagnitude());
                
                if (distanceFromLine > m_distanceFromLine)
                    continue;
                
                extrapolatedCaloHitVector.push_back(pCaloHit);

                if ((pCluster != clusterAssociation.GetUpstreamCluster()) && pCluster != clusterAssociation.GetDownstreamCluster())
                    clusterToCaloHitListMap[pCluster].push_back(pCaloHit);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::IsTrackContinuous(const ClusterAssociation &clusterAssociation, CaloHitVector &extrapolatedCaloHitVector) const
{
    const CartesianVector &upstreamMergePoint(clusterAssociation.GetUpstreamMergePoint()), &downstreamMergePoint(clusterAssociation.GetDownstreamMergePoint());
    const CartesianVector &trackDirection(clusterAssociation.GetConnectingLineDirection());
    const CartesianVector trackStep(trackDirection * m_lineSegmentLength);

    // to handle final segment - merge track remainder with preceding segment and if track remainder is more than half trackstep split into two
    const float trackLength((downstreamMergePoint - upstreamMergePoint).GetMagnitude());
    const unsigned int fullSegments(floor(trackLength / m_lineSegmentLength));
    const float lengthOfTrackRemainder(trackLength - (fullSegments * m_lineSegmentLength));

    // sort hits by projected distance from the upstreamMergePoint
    std::sort(extrapolatedCaloHitVector.begin(), extrapolatedCaloHitVector.end(), SortByDistanceAlongLine(upstreamMergePoint, trackDirection));

    unsigned int segmentsWithoutHits(0);
    CaloHitVector::const_iterator caloHitIter(extrapolatedCaloHitVector.begin());

    for (unsigned int i = 0; i < (fullSegments + 1); ++i)
    {
        // if have run out of hits
        if (caloHitIter == extrapolatedCaloHitVector.end())
        {
            ++segmentsWithoutHits;
            if (segmentsWithoutHits > m_maxTrackGaps)
                return false;
            continue;
        }

        CartesianVector lowerBoundary(upstreamMergePoint + (trackStep * i));
        CartesianVector upperBoundary(upstreamMergePoint + (trackStep * (i + 1)));

        if (i >= fullSegments - 1)
        {
            if (lengthOfTrackRemainder > m_lineSegmentLength * 0.5f)
            {
                lowerBoundary = lowerBoundary - (trackStep * 0.5f * (i - fullSegments + 1.f)) + (trackDirection * 0.5f * (i - fullSegments + 1.f) * lengthOfTrackRemainder);
                upperBoundary = upperBoundary - (trackStep * 0.5f * (i - fullSegments + 2.f)) + (trackDirection * 0.5f * (i - fullSegments + 2.f) * lengthOfTrackRemainder);
            }
            else
            {
                upperBoundary = downstreamMergePoint;
            }
        }

        unsigned int hitsInSegment(0);
        while (this->IsInLineSegment(lowerBoundary, upperBoundary, (*caloHitIter)->GetPositionVector()))
        {
            ++hitsInSegment;
            ++caloHitIter;

            if (caloHitIter == extrapolatedCaloHitVector.end())
                break;
        }

        segmentsWithoutHits = hitsInSegment ? 0 : segmentsWithoutHits + 1;

        if (segmentsWithoutHits > m_maxTrackGaps)
            return false;

        // leave loop early if final merged segment is not split into two
        if (i == (fullSegments - 1) && !(lengthOfTrackRemainder > m_lineSegmentLength * 0.5))
            return true;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::IsInLineSegment(const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary, const CartesianVector &point) const
{
    const float gradient = (-1.0)*(upperBoundary.GetX() - lowerBoundary.GetX()) / (upperBoundary.GetZ() - lowerBoundary.GetZ());
    const float xPointOnUpperLine((point.GetZ() - upperBoundary.GetZ() + gradient*upperBoundary.GetX())/gradient);
    const float xPointOnLowerLine((point.GetZ() - lowerBoundary.GetZ() + gradient*lowerBoundary.GetX())/gradient);

    const CartesianVector upper(xPointOnUpperLine, 0.f, point.GetZ());
    const CartesianVector lower(xPointOnLowerLine, 0.f, point.GetZ());

    if ((point.GetX() > xPointOnUpperLine) && (point.GetX() > xPointOnLowerLine))
        return false;

    if ((point.GetX() < xPointOnUpperLine) && (point.GetX() < xPointOnLowerLine))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::RefineTracks(ClusterAssociation &clusterAssociation, TwoDSlidingFitResultMap &microSlidingFitResultMap,
    TwoDSlidingFitResultMap &macroSlidingFitResultMap, ClusterVector &clusterVector, CaloHitVector &extrapolatedCaloHitVector) const
{
    clusterAssociation.SetUpstreamCluster(this->RefineTrack(clusterAssociation.GetUpstreamCluster(), clusterAssociation.GetUpstreamMergePoint(), extrapolatedCaloHitVector, true,
        microSlidingFitResultMap, macroSlidingFitResultMap, clusterVector));
    clusterAssociation.SetDownstreamCluster(this->RefineTrack(clusterAssociation.GetDownstreamCluster(), clusterAssociation.GetDownstreamMergePoint(), extrapolatedCaloHitVector, false,
        microSlidingFitResultMap, macroSlidingFitResultMap, clusterVector));
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster* TrackInEMShowerAlgorithm::RefineTrack(const Cluster *const pCluster, const CartesianVector &splitPosition, const CaloHitVector &extrapolatedCaloHitVector,
    const bool isUpstream, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap, ClusterVector &clusterVector) const
{
    // Get the split position fitting coordinates
    float rL(0.f), rT(0.f);
    const TwoDSlidingFitResult &microFitResult(microSlidingFitResultMap.at(pCluster));
    microFitResult.GetLocalPosition(splitPosition, rL, rT);

    const TwoDSlidingFitResult macroFitResult(macroSlidingFitResultMap.at(pCluster));
    CartesianVector averageDirection(0.f, 0.f, 0.f);
    macroFitResult.GetGlobalDirection(macroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), averageDirection);

    const float clusterGradient(averageDirection.GetZ()/averageDirection.GetX());           
    const float clusterIntercept(splitPosition.GetZ() - (clusterGradient * splitPosition.GetX()));

    // FRAGMENTATION INITIALISATION
    std::string originalListName, fragmentListName;
    const ClusterList originalClusterList(1, pCluster);    
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

    const Cluster *pMainTrackCluster(nullptr), *pToFragmentCluster(nullptr);
    CaloHitList aboveHitList, belowHitList;
    
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            float thisL(0.f), thisT(0.f);

            microFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), thisL, thisT);

            const bool isAnExtrapolatedHit(std::find(extrapolatedCaloHitVector.begin(), extrapolatedCaloHitVector.end(), pCaloHit) != extrapolatedCaloHitVector.end());
            const bool toRemove(!isAnExtrapolatedHit && (((thisL > rL) && isUpstream) || ((thisL < rL) && !isUpstream)));
                
            if (toRemove)
            {
                CaloHitList &caloHitListToModify((hitPosition.GetZ() > ((clusterGradient * hitPosition.GetX()) + clusterIntercept)) ? aboveHitList : belowHitList);
                caloHitListToModify.push_back(pCaloHit);
            }

            const Cluster *&pClusterToModify(toRemove ? pToFragmentCluster : pMainTrackCluster);
            
            if (pClusterToModify)
            {
                PandoraContentApi::AddToCluster(*this, pClusterToModify, pCaloHit);
            }
            else
            {
                PandoraContentApi::Cluster::Parameters parameters;
                parameters.m_caloHitList.push_back(pCaloHit);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClusterToModify));
            }
        }
    }

    if (!pToFragmentCluster)
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, originalListName, fragmentListName));
        return pCluster;
    }
    else
    {
        this->UpdateForClusterDeletion(pCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
    }

    ClusterList mainTrackClusterList({pMainTrackCluster});
    this->SelectCleanClusters(&mainTrackClusterList, clusterVector);
    
    const bool fragmentAbove(!aboveHitList.empty());
    const bool fragmentBelow(!belowHitList.empty());
    
    if (fragmentAbove || fragmentBelow)
        this->FragmentCluster(pToFragmentCluster, pMainTrackCluster, aboveHitList, belowHitList, fragmentAbove, fragmentBelow);
            
    return pMainTrackCluster;
}

//////////////////////
/*
    std::cout << "NUMBER OF ABOVE HITS: " << aboveHitList.size() << std::endl;
    std::cout << "NUMBER OF BELOW HITS: " << belowHitList.size() << std::endl;
    PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &aboveHitList, "ABOVE TRACK HITS", GREEN);
    PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &belowHitList, "BELOW TRACK HITS", DARKGREEN);

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
*/
//////////////////////

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::FragmentCluster(const Cluster *const pCluster, const Cluster *const pClusterToEnlarge, const CaloHitList &aboveTrackList,
    const CaloHitList &belowTrackList, const bool fragmentAbove, const bool fragmentBelow) const
{
    
    std::vector<std::pair<const CaloHitList*, bool>> fragmentationInfo({{&aboveTrackList, fragmentAbove}, {&belowTrackList, fragmentBelow}});
    
    // FRAGMENTATION INITIALISATION
    std::string originalListName, fragmentListName;
    const ClusterList originalClusterList(1, pCluster);    
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

    ClusterList createdClusters;
    for (auto &fragmentPair : fragmentationInfo)
    {
        if (!fragmentPair.second && !fragmentPair.first->empty())
        {
            const Cluster *pNewCluster(nullptr);
            
            PandoraContentApi::Cluster::Parameters parameters;
            parameters.m_caloHitList = *fragmentPair.first;
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewCluster));
            createdClusters.push_back(pNewCluster);
            continue;
        }
        
        for (const CaloHit *const pCaloHit : *fragmentPair.first) 
        {
            const Cluster *pClosestCluster(nullptr);
            
            if (!createdClusters.empty())
            {
                float closestDistance(std::numeric_limits<float>::max());

                for (const Cluster *const pCreatedCluster : createdClusters)
                {
                    const float separationDistance(LArClusterHelper::GetClosestDistance(pCaloHit->GetPositionVector(), pCreatedCluster));
                    if ((separationDistance < closestDistance) && (separationDistance < m_maxHitDistanceFromCluster))
                    {
                        closestDistance = separationDistance;
                        pClosestCluster = pCreatedCluster;
                    }
                }
            }

            if (pClosestCluster)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClosestCluster, pCaloHit));
            }
            else
            {
                PandoraContentApi::Cluster::Parameters parameters;
                parameters.m_caloHitList.push_back(pCaloHit);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClosestCluster));
                createdClusters.push_back(pClosestCluster);
            }
        }
    }

    // END FRAGMENTATION PROCESS
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName)); 

    for (const Cluster *const pCreatedCluster : createdClusters)
    {
        if (pCreatedCluster->GetNCaloHits() == 1)
            this->AddToNearestCluster(pCreatedCluster, pClusterToEnlarge);
    }

    // ATTN: Do not add clusters that now pass the 'clean cluster' threshold into the cluster vector, since these are unlikely to be a muon track
}

/////////////////////////
/*
            std::cout << "ADDED HIT TO THE NEAREST CLUSTER..." << std::endl;

    for (const Cluster *const pCreatedCluster : createdClusters)
    {
        ClusterList theCreatedCluster({pCreatedCluster});
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &theCreatedCluster, "CREATED CLUSTER", AUTO);
    }
*/
/////////////////////////
//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::MergeHits(const ClusterAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap,
    TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap,  ClusterVector &clusterVector) const
{
    // FIRST DEAL WITH MAIN CLUSTER PAIR
    const Cluster *const pClusterToEnlarge(clusterAssociation.GetUpstreamCluster());
    const Cluster *const pClusterToDelete(clusterAssociation.GetDownstreamCluster());

    this->UpdateForClusterDeletion(pClusterToDelete, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);    
    PandoraContentApi::MergeAndDeleteClusters(*this, pClusterToEnlarge, pClusterToDelete); 

    // NEXT ADD IN ALL OTHER CLUSTERS
    ClusterVector clustersToFragment;
    for (const auto &entry : clusterToCaloHitListMap)
        clustersToFragment.push_back(entry.first);
        
    if (clustersToFragment.empty())
    {
        this->UpdateAfterClusterModification(pClusterToEnlarge, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
        return;
    }

    std::sort(clustersToFragment.begin(), clustersToFragment.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster : clustersToFragment)
    {
        this->UpdateForClusterDeletion(pCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

        const CaloHitList &caloHitsToMerge(clusterToCaloHitListMap.at(pCluster));

        this->MergeCluster(pCluster, pClusterToEnlarge, caloHitsToMerge, clusterAssociation);
    }
        
    this->UpdateAfterClusterModification(pClusterToEnlarge, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::MergeCluster(const Cluster *const pCluster, const Cluster *const pClusterToEnlarge, const CaloHitList &caloHitsToMerge,
    const ClusterAssociation &clusterAssociation) const
{

    if (pCluster->GetNCaloHits() == caloHitsToMerge.size())
    {
        PandoraContentApi::MergeAndDeleteClusters(*this, pClusterToEnlarge, pCluster);
        return;
    }

    CaloHitList aboveTrackHits;
    CaloHitList belowTrackHits;

    const float connectingLineGradient(clusterAssociation.GetConnectingLineDirection().GetZ()/clusterAssociation.GetConnectingLineDirection().GetX());         
    const float connectingLineIntercept(clusterAssociation.GetUpstreamMergePoint().GetZ() - (connectingLineGradient * clusterAssociation.GetUpstreamMergePoint().GetX()));

    const OrderedCaloHitList orderedCaloHitList(pCluster->GetOrderedCaloHitList()); 
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            const bool isAnExtrapolatedHit(std::find(caloHitsToMerge.begin(), caloHitsToMerge.end(), pCaloHit) != caloHitsToMerge.end());
            if (isAnExtrapolatedHit)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, pCluster, pCaloHit));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClusterToEnlarge, pCaloHit));
            }
            else
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                CaloHitList &positionCaloHitList((hitPosition.GetZ() > ((connectingLineGradient * hitPosition.GetX()) + connectingLineIntercept)) ? aboveTrackHits : belowTrackHits);
                
                positionCaloHitList.push_back(pCaloHit);
            }
        }
    }

    if (pCluster->GetNCaloHits() == 1)
    {
        this->AddToNearestCluster(pCluster, pClusterToEnlarge);
        return;
    }

    const bool fragmentAbove(this->IsClusterRemnantDisconnected(aboveTrackHits)), fragmentBelow(this->IsClusterRemnantDisconnected(belowTrackHits));

    this->FragmentCluster(pCluster, pClusterToEnlarge, aboveTrackHits, belowTrackHits, fragmentAbove, fragmentBelow);
}

////////////////////////
/*
    std::cout << "NUMBER OF ABOVE HITS: " << aboveTrackHits.size() << std::endl;
    std::cout << "NUMBER OF BELOW HITS: " << belowTrackHits.size() << std::endl;
    PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &aboveTrackHits, "ABOVE TRACK HITS", GREEN);
    PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &belowTrackHits, "BELOW TRACK HITS", DARKGREEN);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
*/
////////////////////////

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::IsClusterRemnantDisconnected(const CaloHitList &trackHits) const
{
    if (trackHits.size() < 2)
    {
        return false;
    }
        
    CaloHitVector trackHitsVector(trackHits.begin(), trackHits.end());
    std::sort(trackHitsVector.begin(), trackHitsVector.end(), LArClusterHelper::SortHitsByPosition);
    
    float maxSeparationDistance(-std::numeric_limits<float>::max());
    for (CaloHitVector::const_iterator iter = trackHitsVector.begin(); iter != std::prev(trackHitsVector.end(), 1); ++iter)
    {
        const CaloHit *const pFirstCaloHit(*iter), *const pSecondCaloHit(*(std::next(iter, 1)));
        const float separationDistance(std::sqrt(pFirstCaloHit->GetPositionVector().GetDistanceSquared(pSecondCaloHit->GetPositionVector())));

        if (separationDistance > maxSeparationDistance)
            maxSeparationDistance = separationDistance;
    }

    if (maxSeparationDistance > 5.f)
        return true;

    return false;
}



//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::AddToNearestCluster(const CaloHit *const pCaloHitToMerge, const Cluster *const pParentCluster, const Cluster *const pClusterToEnlarge) const
{
    const ClusterList *pClusterList = NULL;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    const CartesianVector &caloHitToMergePosition(pCaloHitToMerge->GetPositionVector());
    
    float closestDistance(std::numeric_limits<float>::max());
    const Cluster *pClosestCluster(nullptr);
    
    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster == pParentCluster)
            continue;
        
        float separationDistance(LArClusterHelper::GetClosestDistance(caloHitToMergePosition, pCluster));

        if (separationDistance < closestDistance)
        {
            if ((pCluster == pClusterToEnlarge) && (separationDistance > 0.75))
                continue;
            
            pClosestCluster = pCluster;
            closestDistance = separationDistance;
        }
    }

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClusterToEnlarge, pCaloHitToMerge));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::AddToNearestCluster(const Cluster *const pClusterToMerge, const Cluster *const pClusterToEnlarge) const
{
    const ClusterList *pClusterList = NULL;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    float closestDistance(std::numeric_limits<float>::max());
    const Cluster *pClosestCluster(nullptr);

    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster == pClusterToMerge)
            continue;
        
        float separationDistance(LArClusterHelper::GetClosestDistance(pClusterToMerge, pCluster));

        if (separationDistance < closestDistance)
        {
            if ((pCluster == pClusterToEnlarge) && (separationDistance > 0.75))
                continue;
            
            pClosestCluster = pCluster;
            closestDistance = separationDistance;
        }
    }

    PandoraContentApi::MergeAndDeleteClusters(*this, pClosestCluster, pClusterToMerge); 
        
}
//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::UpdateAfterClusterModification(const Cluster *const pCluster, ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
    TwoDSlidingFitResultMap &macroSlidingFitResultMap) const
{
    // Remove from SF maps
    std::vector<TwoDSlidingFitResultMap*> slidingFitResultVector({&microSlidingFitResultMap, &macroSlidingFitResultMap});
    this->RemoveClusterFromSlidingFitResultMaps(pCluster, slidingFitResultVector);

    // Remove from CV
    this->RemoveClusterFromClusterVector(pCluster, clusterVector);

    // Check whether cluster should be added back into the CV
    ClusterList alteredClusterList({pCluster});
    this->SelectCleanClusters(&alteredClusterList, clusterVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::UpdateForClusterDeletion(const Cluster *const pCluster, ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
    TwoDSlidingFitResultMap &macroSlidingFitResultMap) const
{
    // Remove from SF maps
    std::vector<TwoDSlidingFitResultMap*> slidingFitResultVector({&microSlidingFitResultMap, &macroSlidingFitResultMap});
    this->RemoveClusterFromSlidingFitResultMaps(pCluster, slidingFitResultVector);

    // Remove from CV
    this->RemoveClusterFromClusterVector(pCluster, clusterVector);
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
void TrackInEMShowerAlgorithm::RemoveClusterFromSlidingFitResultMaps(const Cluster *const pCluster, std::vector<TwoDSlidingFitResultMap*> &slidingFitResultMapVector) const
{
    for (TwoDSlidingFitResultMap *const pSlidingFitResultMap : slidingFitResultMapVector)
    {
        const TwoDSlidingFitResultMap::const_iterator fitToDelete(pSlidingFitResultMap->find(pCluster));
        if (fitToDelete != pSlidingFitResultMap->end())
            pSlidingFitResultMap->erase(fitToDelete);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::RemoveClusterFromClusterVector(const Cluster *const pCluster, ClusterVector &clusterVector) const
{
    ClusterVector::const_iterator clusterToDelete(std::find(clusterVector.begin(), clusterVector.end(), pCluster));
    if (clusterToDelete != clusterVector.end())
        clusterVector.erase(clusterToDelete);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::UpdateSlidingFitResultMap(const ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
    TwoDSlidingFitResultMap &macroSlidingFitResultMap) const
{
    const float slidingFitPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    
    for (const Cluster *const pCluster : clusterVector)
    {
        const TwoDSlidingFitResultMap::const_iterator microFitResultIter(microSlidingFitResultMap.find(pCluster));
        
        if(microFitResultIter == microSlidingFitResultMap.end())
        {
            try
            {
                (void) microSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch)));
            }
            catch (StatusCodeException &) {}
        }

        const TwoDSlidingFitResultMap::const_iterator macroFitResultIter(macroSlidingFitResultMap.find(pCluster));
        
        if(macroFitResultIter == macroSlidingFitResultMap.end())
        {
            try
            {
                (void) macroSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, TwoDSlidingFitResult(pCluster, 1000000, slidingFitPitch)));
            }
            catch (StatusCodeException &) {}
        }
    }
}


//------------------------------------------------------------------------------------------------------------------------------------------
    
StatusCode TrackInEMShowerAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCaloHits", m_minCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SlidingFitWindow", m_slidingFitWindow));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxPredictedMergePointOffset", m_maxPredictedMergePointOffset));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinSeparationDistance", m_minSeparationDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MergePointMinCosAngleDeviation", m_mergePointMinCosAngleDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterLengthSum", m_minClusterLengthSum));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinDirectionDeviationCosAngle", m_minDirectionDeviationCosAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "DistanceFromLine", m_distanceFromLine));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxTrackGaps", m_maxTrackGaps));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LineSegmentLength", m_lineSegmentLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxHitsToRecluster", m_maxHitsToRecluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinModificationFractionToRecluster", m_minModificationFractionToRecluster));
    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content


    
