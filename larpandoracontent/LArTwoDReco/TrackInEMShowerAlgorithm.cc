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
    
TrackInEMShowerAlgorithm::TrackInEMShowerAlgorithm() :
    m_minCaloHits(25), 
    m_minSeparationDistance(27.f),
    m_maxPredictedMergePointOffset(5.f),
    m_slidingFitWindow(20),
    m_mergePointMinCosAngleDeviation(0.9995f),
    m_minClusterLengthSum(75.f),
    m_minDirectionDeviationCosAngle(0.99f),
    m_distanceFromLine(0.35f), 
    m_maxTrackGaps(2),
    m_lineSegmentLength(7.f),
    m_maxHitDistanceFromCluster(4.f),
    m_maxHitSeparationForConnectedCluster(4.f),
    m_maxDistanceFromMainTrack(0.75f),
    m_maxMainLoopIterations(5),
    m_macroSlidingFitWindow(1000),
    m_stableRegionClusterFraction(0.1f),
    m_useLayer(false)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode TrackInEMShowerAlgorithm::Run()
{
    const ClusterList *pClusterList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pClusterList));

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentList(*this, pCaloHitList));
    
    ClusterVector clusterVector;
    this->SelectCleanClusters(pClusterList, clusterVector);

    TwoDSlidingFitResultMap microSlidingFitResultMap;
    TwoDSlidingFitResultMap macroSlidingFitResultMap;
    this->InitialiseSlidingFitResultMaps(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);

    bool mergeMade(false);
    unsigned int loopIterations(0);
    do
    {
        ++loopIterations;
        
        ClusterAssociation clusterAssociation;
        if(!this->FindBestClusterAssociation(clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap, clusterAssociation))
            break;

        ClusterToCaloHitListMap clusterToCaloHitListMap;
        this->GetExtrapolatedCaloHits(clusterAssociation, pClusterList, clusterToCaloHitListMap);

        if (!this->IsTrackContinuous(clusterAssociation, clusterToCaloHitListMap))
            break;
        
        this->CreateMainTrack(clusterAssociation, clusterToCaloHitListMap, pClusterList, microSlidingFitResultMap, macroSlidingFitResultMap, clusterVector);
        
        mergeMade = true;
        
    }
    while(mergeMade && (loopIterations < m_maxMainLoopIterations));

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
    
    for (const Cluster *const pCluster : clusterVector)
    {
        try
        {
            TwoDSlidingFitResult microSlidingFitResult(pCluster, m_slidingFitWindow, slidingFitPitch);
            TwoDSlidingFitResult macroSlidingFitResult(pCluster, m_macroSlidingFitWindow, slidingFitPitch);

            (void) microSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, microSlidingFitResult));
            (void) macroSlidingFitResultMap.insert(TwoDSlidingFitResultMap::value_type(pCluster, macroSlidingFitResult));
        }
        catch (StatusCodeException &) {}
    }
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::FindBestClusterAssociation(const ClusterVector &clusterVector, const TwoDSlidingFitResultMap &microSlidingFitResultMap,
    const TwoDSlidingFitResultMap &macroSlidingFitResultMap, ClusterAssociation &clusterAssociation) const
{
    bool foundAssociation(false);
    float maxLength(0.f);

    // ATTN: Find the associated cluster pair with the longest length sum
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

            // ATTN: Ensure that clusters are not contained within one another
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

            if (isCurrentUpstream)
            {
                clusterAssociation = ClusterAssociation(pCurrentCluster, pTestCluster, currentMergePoint, currentMergeDirection, testMergePoint, testMergeDirection);
            }
            else
            {
                clusterAssociation = ClusterAssociation(pTestCluster, pCurrentCluster, testMergePoint, testMergeDirection, currentMergePoint, currentMergeDirection);
            }

            foundAssociation = true;
            maxLength = lengthSum;
        }
    }

    return foundAssociation;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::GetClusterMergingCoordinates(const TwoDSlidingFitResult &clusterMicroFitResult, const TwoDSlidingFitResult &clusterMacroFitResult,
    const TwoDSlidingFitResult &associatedMacroFitResult, CartesianVector &clusterMergePosition, CartesianVector &clusterMergeDirection, const bool isUpstream) const
{
    CartesianVector clusterAverageDirection(0.f, 0.f, 0.f);
    CartesianVector associatedAverageDirection(0.f, 0.f, 0.f);
    clusterMacroFitResult.GetGlobalDirection(clusterMacroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), clusterAverageDirection);
    associatedMacroFitResult.GetGlobalDirection(associatedMacroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), associatedAverageDirection);    
    
    const LayerFitResultMap &clusterMicroLayerFitResultMap(clusterMicroFitResult.GetLayerFitResultMap());
    const int startLayer(isUpstream ? clusterMicroFitResult.GetMaxLayer() : clusterMicroFitResult.GetMinLayer());
    const int endLayer(isUpstream ? clusterMicroFitResult.GetMinLayer() : clusterMicroFitResult.GetMaxLayer());
    const int loopTerminationLayer(endLayer + (isUpstream ? -1 : 1));
    const int step(isUpstream ? -1 : 1);

    // ATTN: Search for stable region for which the local layer gradient agrees well with associated cluster global gradient
    unsigned int gradientStabilityWindow;
    if (m_useLayer)
    {
        gradientStabilityWindow = std::ceil(clusterMicroLayerFitResultMap.size() *  m_stableRegionClusterFraction);
    }
    else
    {
        gradientStabilityWindow = std::ceil(clusterMicroFitResult.GetCluster()->GetNCaloHits() * m_stableRegionClusterFraction);
    }
    unsigned int goodLayerCount(0);      

    for (int i = startLayer; i != loopTerminationLayer; i += step)
    {
        const auto microIter(clusterMicroLayerFitResultMap.find(i));

        if (microIter == clusterMicroLayerFitResultMap.end())
            continue;
        
        CartesianVector microDirection(0.f, 0.f, 0.f);
        clusterMicroFitResult.GetGlobalDirection(microIter->second.GetGradient(), microDirection);
    
        const float cosDirectionOpeningAngle(microDirection.GetCosOpeningAngle(associatedAverageDirection));
        if (cosDirectionOpeningAngle > m_mergePointMinCosAngleDeviation)
        {
            // ATTN: Set merge point and direction as that at the first layer in the stable region
            if (goodLayerCount == 0)
            {
                // ATTN: Cluster direction vectors must point to one another
                clusterMergeDirection = clusterAverageDirection * (isUpstream ? 1.f : -1.f);
                clusterMicroFitResult.GetGlobalFitPosition(microIter->second.GetL(), clusterMergePosition);
            }
            
            ++goodLayerCount;
        }
        else
        {
            goodLayerCount = 0;
        }

        if (goodLayerCount > gradientStabilityWindow)
            break;
        
        // ATTN: Abort merging process have not found a stable region 
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

    const float separationDistance(std::sqrt(upstreamPoint.GetDistanceSquared(downstreamPoint)));
    if (separationDistance < m_minSeparationDistance)
        return false;
    
    // ATTN: Cluster directions point towards one another
    if (upstreamDirection.GetCosOpeningAngle(downstreamDirection * (-1.f)) < m_minDirectionDeviationCosAngle)
        return false;
    
    const CartesianVector predictedDownstreamPoint(upstreamPoint + (upstreamDirection * separationDistance));
    const float predictedDownstreamOffset((predictedDownstreamPoint - downstreamPoint).GetMagnitude());
    
    if (predictedDownstreamOffset > m_maxPredictedMergePointOffset)
        return false;

    const CartesianVector predictedUpstreamPoint(downstreamPoint + (downstreamDirection * separationDistance));
    const float predictedUpstreamOffset((predictedUpstreamPoint - upstreamPoint).GetMagnitude());

    if (predictedUpstreamOffset > m_maxPredictedMergePointOffset)
        return false;
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::GetExtrapolatedCaloHits(const ClusterAssociation &clusterAssociation, const ClusterList *const pClusterList,
    ClusterToCaloHitListMap &clusterToCaloHitListMap) const
{
    const CartesianVector &upstreamPoint(clusterAssociation.GetUpstreamMergePoint()), &downstreamPoint(clusterAssociation.GetDownstreamMergePoint());
    const float minX(std::min(upstreamPoint.GetX(), downstreamPoint.GetX())), maxX(std::max(upstreamPoint.GetX(), downstreamPoint.GetX()));

    for (const Cluster *const pCluster : *pClusterList) 
    {
        const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
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
                
                clusterToCaloHitListMap[pCluster].push_back(pCaloHit);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::IsTrackContinuous(const ClusterAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap) const
{
    if (m_lineSegmentLength < std::numeric_limits<float>::epsilon())
    {
        std::cout << "TrackInEMShowerAlgorithm: Line segment length must be positive and nonzero" << std::endl;
        return false;
    }

    const CartesianVector &upstreamMergePoint(clusterAssociation.GetUpstreamMergePoint()), &downstreamMergePoint(clusterAssociation.GetDownstreamMergePoint());
    const CartesianVector &trackDirection(clusterAssociation.GetConnectingLineDirection());
    const CartesianVector trackStep(trackDirection * m_lineSegmentLength);

    // To handle final segment - merge track remainder with preceding segment and if track remainder was more than half trackstep split into two
    const float trackLength((downstreamMergePoint - upstreamMergePoint).GetMagnitude());
    const int fullSegments(std::floor(trackLength / m_lineSegmentLength));
    const float lengthOfTrackRemainder(trackLength - (fullSegments * m_lineSegmentLength));

    // Collect extrapolated calo hits and sort by projected distance from the upstreamMergePoint
    CaloHitVector extrapolatedCaloHitVector;
    for (const auto &entry : clusterToCaloHitListMap)
        extrapolatedCaloHitVector.insert(extrapolatedCaloHitVector.begin(), entry.second.begin(), entry.second.end());
    
    std::sort(extrapolatedCaloHitVector.begin(), extrapolatedCaloHitVector.end(), SortByDistanceAlongLine(upstreamMergePoint, trackDirection));

    unsigned int segmentsWithoutHits(0);
    CaloHitVector::const_iterator caloHitIter(extrapolatedCaloHitVector.begin());
    
    for (int i = 0; i < (fullSegments + 1); ++i)
    {
        // if have run out of hits
        if (caloHitIter == extrapolatedCaloHitVector.end())
        {
            ++segmentsWithoutHits;
            
            if (segmentsWithoutHits > m_maxTrackGaps)
                return false;
            
            continue;
        }

        CartesianVector lowerBoundary(upstreamMergePoint + (trackStep * static_cast<float>(i)));
        CartesianVector upperBoundary(upstreamMergePoint + (trackStep * static_cast<float>(i + 1)));

        if (i >= fullSegments - 1)
        {
            if (lengthOfTrackRemainder > m_lineSegmentLength * 0.5f)
            {
                lowerBoundary = lowerBoundary - (trackStep * 0.5f * (static_cast<float>(i - fullSegments) + 1.f)) + (trackDirection * 0.5f * (static_cast<float>(i - fullSegments) + 1.f)
                    * lengthOfTrackRemainder);
                upperBoundary = upperBoundary - (trackStep * 0.5f * (static_cast<float>(i - fullSegments) + 2.f)) + (trackDirection * 0.5f * (static_cast<float>(i - fullSegments) + 2.f)
                    * lengthOfTrackRemainder);
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
        if (i == (fullSegments - 1) && !(lengthOfTrackRemainder > m_lineSegmentLength * 0.5f))
            return true;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::IsInLineSegment(const CartesianVector &lowerBoundary, const CartesianVector &upperBoundary, const CartesianVector &point) const
{
    const float segmentBoundaryGradient = (-1.f) * (upperBoundary.GetX() - lowerBoundary.GetX()) / (upperBoundary.GetZ() - lowerBoundary.GetZ());
    const float xPointOnUpperLine((point.GetZ() - upperBoundary.GetZ()) / segmentBoundaryGradient + upperBoundary.GetX());
    const float xPointOnLowerLine((point.GetZ() - lowerBoundary.GetZ()) / segmentBoundaryGradient + lowerBoundary.GetX());
    const CartesianVector upper(xPointOnUpperLine, 0.f, point.GetZ());
    const CartesianVector lower(xPointOnLowerLine, 0.f, point.GetZ());

    if ((point.GetX() > xPointOnUpperLine) && (point.GetX() > xPointOnLowerLine))
        return false;

    if ((point.GetX() < xPointOnUpperLine) && (point.GetX() < xPointOnLowerLine))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::CreateMainTrack(const ClusterAssociation &clusterAssociation, const ClusterToCaloHitListMap &clusterToCaloHitListMap, const ClusterList *const pClusterList, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap, ClusterVector &clusterVector) const
{
    // Determine the shower clusters which contain hits that belong to the main track
    ClusterVector showerClustersToFragment;
    for (const auto &entry : clusterToCaloHitListMap)
    {
        if ((entry.first != clusterAssociation.GetUpstreamCluster()) && (entry.first != clusterAssociation.GetDownstreamCluster()))
            showerClustersToFragment.push_back(entry.first);
    }

    std::sort(showerClustersToFragment.begin(), showerClustersToFragment.end(), LArClusterHelper::SortByNHits);

    ClusterList remnantClusterList;
    
    const Cluster *const pMainTrackCluster(this->RemoveOffAxisHitsFromTrack(clusterAssociation.GetUpstreamCluster(), clusterAssociation.GetUpstreamMergePoint(),
        true, clusterToCaloHitListMap, remnantClusterList, microSlidingFitResultMap, macroSlidingFitResultMap, clusterVector));
    const Cluster *const pClusterToDelete(this->RemoveOffAxisHitsFromTrack(clusterAssociation.GetDownstreamCluster(), clusterAssociation.GetDownstreamMergePoint(),
        false, clusterToCaloHitListMap, remnantClusterList, microSlidingFitResultMap, macroSlidingFitResultMap, clusterVector));

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pMainTrackCluster, pClusterToDelete));

    for (const Cluster *const pShowerCluster : showerClustersToFragment)
    {
        const CaloHitList &caloHitsToMerge(clusterToCaloHitListMap.at(pShowerCluster));
        
        this->UpdateForClusterDeletion(pShowerCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
        this->AddHitsToMainTrack(pShowerCluster, pMainTrackCluster, caloHitsToMerge, clusterAssociation, remnantClusterList);
    }

    this->ProcessRemnantClusters(remnantClusterList, pMainTrackCluster, pClusterList);
    this->UpdateAfterMainTrackCreation(pMainTrackCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *TrackInEMShowerAlgorithm::RemoveOffAxisHitsFromTrack(const Cluster *const pCluster, const CartesianVector &splitPosition, const bool isUpstream,
    const ClusterToCaloHitListMap &clusterToCaloHitListMap, ClusterList &remnantClusterList, TwoDSlidingFitResultMap &microSlidingFitResultMap, TwoDSlidingFitResultMap &macroSlidingFitResultMap,
    ClusterVector &clusterVector) const
{
    const auto extrapolatedCaloHitIter(clusterToCaloHitListMap.find(pCluster));

    if (extrapolatedCaloHitIter == clusterToCaloHitListMap.end())
    {
        this->UpdateForClusterDeletion(pCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
        return pCluster;
    }
    
    float rL(0.f), rT(0.f);
    const TwoDSlidingFitResult &microFitResult(microSlidingFitResultMap.at(pCluster));
    microFitResult.GetLocalPosition(splitPosition, rL, rT);

    const TwoDSlidingFitResult macroFitResult(macroSlidingFitResultMap.at(pCluster));
    CartesianVector averageDirection(0.f, 0.f, 0.f);
    macroFitResult.GetGlobalDirection(macroFitResult.GetLayerFitResultMap().begin()->second.GetGradient(), averageDirection);

    const bool isVertical(std::fabs(averageDirection.GetX()) < std::numeric_limits<float>::epsilon());
    const float clusterGradient(isVertical ? 0.f : averageDirection.GetZ() / averageDirection.GetX());           
    const float clusterIntercept(isVertical ? splitPosition.GetX() : splitPosition.GetZ() - (clusterGradient * splitPosition.GetX()));

    // Fragmentation initialisation
    std::string originalListName, fragmentListName;
    const ClusterList originalClusterList(1, pCluster);    
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

    const Cluster *pMainTrackCluster(nullptr), *pAboveCluster(nullptr), *pBelowCluster(nullptr);

    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            float thisL(0.f), thisT(0.f);

            microFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), thisL, thisT);

            const bool isAnExtrapolatedHit(std::find(extrapolatedCaloHitIter->second.begin(), extrapolatedCaloHitIter->second.end(), pCaloHit) != extrapolatedCaloHitIter->second.end());
            const bool isAbove(((clusterGradient * hitPosition.GetX()) + clusterIntercept) < (isVertical ? hitPosition.GetX() : hitPosition.GetZ()));
            const bool isToRemove(!isAnExtrapolatedHit && (((thisL > rL) && isUpstream) || ((thisL < rL) && !isUpstream)));
                
            const Cluster *&pClusterToModify(isToRemove ? (isAbove ? pAboveCluster : pBelowCluster) : pMainTrackCluster);
            
            if (pClusterToModify)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClusterToModify, pCaloHit));
            }
            else
            {
                PandoraContentApi::Cluster::Parameters parameters;
                parameters.m_caloHitList.push_back(pCaloHit);
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClusterToModify));

                if (pClusterToModify != pMainTrackCluster)
                    remnantClusterList.push_back(pClusterToModify);  
            }
        }
    }

    this->UpdateForClusterDeletion(pCluster, clusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
    
    // End fragmentation
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
    return pMainTrackCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::AddHitsToMainTrack(const Cluster *const pShowerCluster, const Cluster *const pMainTrackCluster, const CaloHitList &caloHitsToMerge,
    const ClusterAssociation &clusterAssociation, ClusterList &remnantClusterList) const
{
    if (pShowerCluster->GetNCaloHits() == caloHitsToMerge.size())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pMainTrackCluster, pShowerCluster));
        return;
    }

    // Fragmentation initialisation
    std::string originalListName, fragmentListName;
    const ClusterList originalClusterList(1, pShowerCluster);    
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

    const Cluster *pAboveCluster(nullptr), *pBelowCluster(nullptr);
    
    const float connectingLineGradient(clusterAssociation.GetConnectingLineDirection().GetZ()/clusterAssociation.GetConnectingLineDirection().GetX());         
    const float connectingLineIntercept(clusterAssociation.GetUpstreamMergePoint().GetZ() - (connectingLineGradient * clusterAssociation.GetUpstreamMergePoint().GetX()));

    const OrderedCaloHitList orderedCaloHitList(pShowerCluster->GetOrderedCaloHitList()); 
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            const bool isAnExtrapolatedHit(std::find(caloHitsToMerge.begin(), caloHitsToMerge.end(), pCaloHit) != caloHitsToMerge.end());
            if (isAnExtrapolatedHit)
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, pShowerCluster, pCaloHit));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pMainTrackCluster, pCaloHit));
            }
            else
            {
                const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
                const bool isAbove(hitPosition.GetZ() > ((connectingLineGradient * hitPosition.GetX()) + connectingLineIntercept));
                const Cluster *&pClusterToModify(isAbove ? pAboveCluster : pBelowCluster);
            
                if (pClusterToModify)
                {
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pClusterToModify, pCaloHit));
                }
                else
                {
                    PandoraContentApi::Cluster::Parameters parameters;
                    parameters.m_caloHitList.push_back(pCaloHit);
                    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pClusterToModify));

                    remnantClusterList.push_back(pClusterToModify);
                }
            }
        }
    }

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool TrackInEMShowerAlgorithm::IsClusterRemnantDisconnected(const Cluster *const pRemnantCluster) const
{
    if (pRemnantCluster->GetNCaloHits() == 1)
        return false;

    unsigned int count(0);
    float maxSeparationDistance(-std::numeric_limits<float>::max());

    const OrderedCaloHitList &orderedCaloHitList(pRemnantCluster->GetOrderedCaloHitList());
    CartesianVector previousHitPosition((*orderedCaloHitList.begin()->second->begin())->GetPositionVector());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
        {
            ++count;
            
            if (count == 1)
                continue;
            
            const float separationDistance(std::sqrt(pCaloHit->GetPositionVector().GetDistanceSquared(previousHitPosition)));

            if (separationDistance > maxSeparationDistance)
                maxSeparationDistance = separationDistance;

            previousHitPosition = pCaloHit->GetPositionVector();
        }
    }

    if (maxSeparationDistance > m_maxHitSeparationForConnectedCluster)
        return true;
    
    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::FragmentRemnantCluster(const Cluster *const pRemnantCluster, ClusterList &fragmentedClusterList) const
{
    // Fragmentation initialisation
    std::string originalListName, fragmentListName;
    const ClusterList originalClusterList(1, pRemnantCluster);    
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::InitializeFragmentation(*this, originalClusterList, originalListName, fragmentListName));

    ClusterList createdClusters;

    const OrderedCaloHitList &orderedCaloHitList(pRemnantCluster->GetOrderedCaloHitList());
    for (const OrderedCaloHitList::value_type &mapEntry : orderedCaloHitList)
    {
        for (const CaloHit *const pCaloHit : *mapEntry.second) 
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

    // End fragmentation
    if (createdClusters.empty())
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, originalListName, fragmentListName));
        fragmentedClusterList.push_back(pRemnantCluster);
    }
    else
    {
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, fragmentListName, originalListName));
        fragmentedClusterList.insert(fragmentedClusterList.begin(), createdClusters.begin(), createdClusters.end());
    }


}


//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::AddToNearestCluster(const Cluster *const pClusterToMerge, const Cluster *const pClusterToEnlarge, const ClusterList *const pClusterList) const
{
    float closestDistance(std::numeric_limits<float>::max());
    const Cluster *pClosestCluster(nullptr);

    for (const Cluster *const pCluster : *pClusterList)
    {
        if (pCluster == pClusterToMerge)
            continue;
        
        float separationDistance(LArClusterHelper::GetClosestDistance(pClusterToMerge, pCluster));

        if (separationDistance < closestDistance)
        {
            if ((pCluster == pClusterToEnlarge) && (separationDistance > m_maxDistanceFromMainTrack))
                continue;
            
            pClosestCluster = pCluster;
            closestDistance = separationDistance;
        }
    }

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::MergeAndDeleteClusters(*this, pClosestCluster, pClusterToMerge));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::ProcessRemnantClusters(const ClusterList &remnantClusterList, const Cluster *const pMainTrackCluster, const ClusterList *const pClusterList) const
{
    for (const Cluster *const pRemnantCluster : remnantClusterList)
    {
        if (pRemnantCluster->GetNCaloHits() == 1)
        {
            this->AddToNearestCluster(pRemnantCluster, pMainTrackCluster, pClusterList);
            continue;
        }
        
        if (this->IsClusterRemnantDisconnected(pRemnantCluster))
        {
            ClusterList fragmentedClusters;
            this->FragmentRemnantCluster(pRemnantCluster, fragmentedClusters);
            
            for (const Cluster *const pFragmentedCluster : fragmentedClusters)
            {
                if (pFragmentedCluster->GetNCaloHits() == 1)
                    this->AddToNearestCluster(pFragmentedCluster, pMainTrackCluster, pClusterList);
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::UpdateForClusterDeletion(const Cluster *const pCluster, ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
    TwoDSlidingFitResultMap &macroSlidingFitResultMap) const
{
    SlidingFitResultMapVector slidingFitResultMapVector({&microSlidingFitResultMap, &macroSlidingFitResultMap});
    
    for (TwoDSlidingFitResultMap *const pSlidingFitResultMap : slidingFitResultMapVector)
    {
        const TwoDSlidingFitResultMap::const_iterator fitToDelete(pSlidingFitResultMap->find(pCluster));
        if (fitToDelete != pSlidingFitResultMap->end())
            pSlidingFitResultMap->erase(fitToDelete);
    }

    ClusterVector::const_iterator clusterToDelete(std::find(clusterVector.begin(), clusterVector.end(), pCluster));
    if (clusterToDelete != clusterVector.end())
        clusterVector.erase(clusterToDelete);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void TrackInEMShowerAlgorithm::UpdateAfterMainTrackCreation(const Cluster *const pMainTrackCluster, ClusterVector &clusterVector, TwoDSlidingFitResultMap &microSlidingFitResultMap,
    TwoDSlidingFitResultMap &macroSlidingFitResultMap) const
{
    clusterVector.push_back(pMainTrackCluster);

    const ClusterVector mainTrackClusterVector({pMainTrackCluster});
    this->InitialiseSlidingFitResultMaps(mainTrackClusterVector, microSlidingFitResultMap, macroSlidingFitResultMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------
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
        "MaxHitSeparationForConnectedCluster", m_maxHitSeparationForConnectedCluster));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxDistanceFromMainTrack", m_maxDistanceFromMainTrack));        

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxMainLoopIterations", m_maxMainLoopIterations));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MacroSlidingFitWindow", m_macroSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "StableRegionClusterFraction", m_stableRegionClusterFraction));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "UseLayer", m_useLayer));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
