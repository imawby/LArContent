/**
 *  @file   LArContent/src/LArHelpers/LArClusterHelper.cc
 * 
 *  @brief  Implementation of the cluster helper class.
 * 
 *  $Log: $
 */

#include "Helpers/ClusterHelper.h"
#include "Helpers/XmlHelper.h"

#include "LArHelpers/LArClusterHelper.h"

#include "LArPseudoLayerCalculator.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar
{

void LArClusterHelper::LArTwoDSlidingFit(const pandora::Cluster *const pCluster, const unsigned int layerFitHalfWindow, TwoDSlidingFitResult &twoDSlidingFitResult)
{
    ClusterHelper::ClusterFitResult clusterFitResult;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, ClusterHelper::FitFullCluster(pCluster, clusterFitResult));

    const CartesianVector axisDirection(clusterFitResult.GetDirection());
    const CartesianVector innerCentroid(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()));
    const CartesianVector axisIntercept(clusterFitResult.GetIntercept() + axisDirection * (axisDirection.GetDotProduct(innerCentroid - clusterFitResult.GetIntercept())));

    LArClusterHelper::LArTwoDSlidingFit(pCluster, layerFitHalfWindow, axisIntercept, axisDirection, twoDSlidingFitResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::LArTwoDSlidingXZFit(const pandora::Cluster *const pCluster, const unsigned int layerFitHalfWindow, TwoDSlidingFitResult &twoDSlidingFitResult)
{
    const CartesianVector axisIntercept(0.f, 0.f, 0.f);
    const CartesianVector axisDirection(0.f, 0.f, 1.f);

    LArClusterHelper::LArTwoDSlidingFit(pCluster, layerFitHalfWindow, axisIntercept, axisDirection, twoDSlidingFitResult);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::LArTwoDSlidingFit(const Cluster *const pCluster, const unsigned int layerFitHalfWindow, const CartesianVector &axisIntercept,
    const CartesianVector &axisDirection, TwoDSlidingFitResult &twoDSlidingFitResult)
{
    if ((std::fabs(axisIntercept.GetY()) > std::numeric_limits<float>::epsilon()) || (std::fabs(axisDirection.GetY()) > std::numeric_limits<float>::epsilon()))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    twoDSlidingFitResult.m_pCluster = pCluster;
    twoDSlidingFitResult.m_layerFitHalfWindow = layerFitHalfWindow;
    twoDSlidingFitResult.m_axisIntercept = axisIntercept;
    twoDSlidingFitResult.m_axisDirection = axisDirection;
    TwoDSlidingFitResult::LayerFitResultMap &layerFitResultMap(twoDSlidingFitResult.m_layerFitResultMap);
    TwoDSlidingFitResult::LayerFitContributionMap &layerFitContributionMap(twoDSlidingFitResult.m_layerFitContributionMap);

    // Identify fit contributions
    int innerLayer(std::numeric_limits<int>::max()), outerLayer(0);
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            float rL(0.f), rT(0.f);
            twoDSlidingFitResult.GetLocalCoordinates((*hitIter)->GetPositionVector(), rL, rT);
            const int layer(twoDSlidingFitResult.GetLayer(rL));

            layerFitContributionMap[layer].AddPoint(rL, rT);

            if (layer < innerLayer)
                innerLayer = layer;

            if (layer > outerLayer)
                outerLayer = layer;
        }
    }

    // Summation for first layers
    unsigned int slidingNPoints(0);
    double slidingSumT(0.), slidingSumL(0.), slidingSumTT(0.), slidingSumLT(0.), slidingSumLL(0.);

    for (int iLayer = innerLayer; iLayer < static_cast<int>(innerLayer + layerFitHalfWindow); ++iLayer)
    {
        TwoDSlidingFitResult::LayerFitContributionMap::const_iterator lyrIter = layerFitContributionMap.find(iLayer);

        if (layerFitContributionMap.end() != lyrIter)
        {
            slidingSumT += lyrIter->second.GetSumT();
            slidingSumL += lyrIter->second.GetSumL();
            slidingSumTT += lyrIter->second.GetSumTT();
            slidingSumLT += lyrIter->second.GetSumLT();
            slidingSumLL += lyrIter->second.GetSumLL();
            slidingNPoints += lyrIter->second.GetNPoints();
        }
    }

    // Sliding fit
    for (int iLayer = innerLayer; iLayer <= outerLayer; ++iLayer)
    {

        const int fwdLayer(iLayer + layerFitHalfWindow);
        TwoDSlidingFitResult::LayerFitContributionMap::const_iterator fwdIter = layerFitContributionMap.find(fwdLayer);

        if (layerFitContributionMap.end() != fwdIter)
        {
            slidingSumT += fwdIter->second.GetSumT();
            slidingSumL += fwdIter->second.GetSumL();
            slidingSumTT += fwdIter->second.GetSumTT();
            slidingSumLT += fwdIter->second.GetSumLT();
            slidingSumLL += fwdIter->second.GetSumLL();
            slidingNPoints += fwdIter->second.GetNPoints();
        }

        const int bwdLayer(iLayer - layerFitHalfWindow - 1);
        TwoDSlidingFitResult::LayerFitContributionMap::const_iterator bwdIter = layerFitContributionMap.find(bwdLayer);

        if (layerFitContributionMap.end() != bwdIter)
        {
            slidingSumT -= bwdIter->second.GetSumT();
            slidingSumL -= bwdIter->second.GetSumL();
            slidingSumTT -= bwdIter->second.GetSumTT();
            slidingSumLT -= bwdIter->second.GetSumLT();
            slidingSumLL -= bwdIter->second.GetSumLL();
            slidingNPoints -= bwdIter->second.GetNPoints();
        }

        if (slidingNPoints > 0)
        {
            const double denominator(slidingSumLL - slidingSumL * slidingSumL / static_cast<double>(slidingNPoints));

            if (std::fabs(denominator) < std::numeric_limits<double>::epsilon())
                continue;

            const double gradient((slidingSumLT - slidingSumL * slidingSumT / static_cast<double>(slidingNPoints)) / denominator);
            const double intercept((slidingSumLL * slidingSumT / static_cast<double>(slidingNPoints) - slidingSumL * slidingSumLT / static_cast<double>(slidingNPoints)) / denominator);

            const double l(twoDSlidingFitResult.GetL(iLayer));
            const double fitT(intercept + gradient * l);

            const double variance((slidingSumTT - 2. * intercept * slidingSumT - 2. * gradient * slidingSumLT + intercept * intercept * static_cast<double>(slidingNPoints) + 2. * gradient * intercept * slidingSumL + gradient * gradient * slidingSumLL) / (1. + gradient * gradient));
            const double rms(std::sqrt(variance / static_cast<double>(slidingNPoints)));

            const TwoDSlidingFitResult::TwoDSlidingFitResult::LayerFitResult layerFitResult(l, fitT, gradient, rms);
            (void) layerFitResultMap.insert(TwoDSlidingFitResult::LayerFitResultMap::value_type(iLayer, layerFitResult));
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::LArTrackWidth(const Cluster *const pCluster)
{
    TwoDSlidingFitResult twoDSlidingFitResult;
    LArClusterHelper::LArTwoDSlidingFit(pCluster, m_layerFitHalfWindow, twoDSlidingFitResult);

    return twoDSlidingFitResult.GetSlidingFitWidth();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetLengthSquared(const Cluster* const pCluster)
{
    const CartesianVector innerCentroid(pCluster->GetCentroid(pCluster->GetInnerPseudoLayer()));
    const CartesianVector outerCentroid(pCluster->GetCentroid(pCluster->GetOuterPseudoLayer()));

    return (outerCentroid - innerCentroid).GetMagnitudeSquared();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetLength(const Cluster* const pCluster)
{
    return std::sqrt(LArClusterHelper::GetLengthSquared(pCluster));
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetEnergyFromLength(const Cluster* const pCluster)
{
    static const float dEdX(0.002f); // approximately 2 MeV/cm

    return (dEdX * LArClusterHelper::GetLength(pCluster));
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArClusterHelper::GetLayerSpan(const Cluster* const pCluster) 
{
    return (1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetLayerOccupancy(const Cluster *const pCluster)
{
    const unsigned int nOccupiedLayers(pCluster->GetOrderedCaloHitList().size());
    const unsigned int nLayers(1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer());

    if (nLayers > 0)
        return (static_cast<float>(nOccupiedLayers) / static_cast<float>(nLayers));

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetLayerOccupancy(const Cluster *const pCluster1, const Cluster *const pCluster2)
{
    const unsigned int nOccupiedLayers(pCluster1->GetOrderedCaloHitList().size() + pCluster2->GetOrderedCaloHitList().size());
    const unsigned int nLayers(1 + std::max(pCluster1->GetOuterPseudoLayer(), pCluster2->GetOuterPseudoLayer()) -
        std::min(pCluster1->GetInnerPseudoLayer(), pCluster2->GetInnerPseudoLayer()));

    if (nLayers > 0)
        return (static_cast<float>(nOccupiedLayers) / static_cast<float>(nLayers));

    return 0.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const Cluster *const pCluster1, const Cluster *const pCluster2)
{
    const OrderedCaloHitList &orderedCaloHitList1(pCluster1->GetOrderedCaloHitList());
    const OrderedCaloHitList &orderedCaloHitList2(pCluster2->GetOrderedCaloHitList());

    bool distanceFound(false);
    float closestDistanceSquared(std::numeric_limits<float>::max());

    for (OrderedCaloHitList::const_iterator iter1 = orderedCaloHitList1.begin(), iterEnd1 = orderedCaloHitList1.end(); iter1 != iterEnd1; ++iter1)
    {
        const CartesianVector position1(pCluster1->GetCentroid(iter1->first));
        
        for (OrderedCaloHitList::const_iterator iter2 = orderedCaloHitList2.begin(), iterEnd2 = orderedCaloHitList2.end(); iter2 != iterEnd2; ++iter2)
        {
            const CartesianVector position2(pCluster2->GetCentroid(iter2->first));
            const float distanceSquared((position2 - position1).GetMagnitudeSquared());

            if (distanceSquared < closestDistanceSquared)
            {
                closestDistanceSquared = distanceSquared;
                distanceFound = true;
            }
        }
    }

    if (distanceFound)
        return std::sqrt(closestDistanceSquared);

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::GetClosestDistance(const CartesianVector &position, const Cluster *const pCluster)
{
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    bool distanceFound(false);
    float closestDistanceSquared(std::numeric_limits<float>::max());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(), iterEnd = orderedCaloHitList.end(); iter != iterEnd; ++iter)
    {
        const CartesianVector layerCentroid(pCluster->GetCentroid(iter->first));
        const float distanceSquared((layerCentroid - position).GetMagnitudeSquared());

        if (distanceSquared < closestDistanceSquared)
        {
            closestDistanceSquared = distanceSquared;
            distanceFound = true;
        }
    }
    
    if (distanceFound)
        return std::sqrt(closestDistanceSquared);

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::GetListOfCleanClusters(const ClusterQuality method, const ClusterList *const pClusterList, ClusterVector &clusterVector)
{
    for (ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter)
    {
        Cluster *pCluster = *iter;

        if (LArClusterHelper::IsCleanCluster(method,pCluster))
            clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

bool LArClusterHelper::IsCleanCluster(const ClusterQuality method, const Cluster* const pCluster)
{
    // TODO: NEED TO RATIONALISE ALL OF THIS SOMEHOW.... HELP!
    if (method == METHOD_A)
        return LArClusterHelper::IsCleanCluster_MethodA(pCluster);

    else if (method == METHOD_B)
        return LArClusterHelper::IsCleanCluster_MethodB(pCluster);
    
    else if (method == METHOD_C)
        return LArClusterHelper::IsCleanCluster_MethodC(pCluster);

    else if (method == METHOD_D)
        return LArClusterHelper::IsCleanCluster_MethodD(pCluster);

    else throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------ 

bool LArClusterHelper::IsCleanCluster_MethodA(const Cluster *const pCluster )
{
    // from ClusterAssociationAlgorithm
    if (LArClusterHelper::GetLayerSpan(pCluster) < 4)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::IsCleanCluster_MethodB(const Cluster *const pCluster)
{
    // from ClusterExtensionAlgorithm, ClusterMergingAlgorithm, VertexFindingAlgorithm
    if (LArClusterHelper::GetLayerSpan(pCluster) < 15)
        return false;

    if (LArClusterHelper::GetLengthSquared(pCluster) < 25.f)
        return false;

    if (LArClusterHelper::GetLayerOccupancy(pCluster) < 0.75f)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::IsCleanCluster_MethodC(const Cluster *const pCluster)
{
    // from VertexSeedFindingAlgorithm
    if (LArClusterHelper::GetLayerSpan(pCluster) < 10)
        return false;

    if (LArClusterHelper::GetLengthSquared(pCluster) < 1.f)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::IsCleanCluster_MethodD(const Cluster *const pCluster)
{
    // from SeedGrowingAlgorithm
    if (pCluster->GetNCaloHits() < 10)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByInnerLayer(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int innerLayerLhs(pLhs->GetInnerPseudoLayer());
    const unsigned int innerLayerRhs(pRhs->GetInnerPseudoLayer());

    if( innerLayerLhs != innerLayerRhs )
      return (innerLayerLhs < innerLayerRhs);

    // Use SortByNOccupiedLayers method to resolve ties
    return SortByNOccupiedLayers(pLhs,pRhs);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByNOccupiedLayers(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int nOccupiedLayersLhs(pLhs->GetOrderedCaloHitList().size());
    const unsigned int nOccupiedLayersRhs(pRhs->GetOrderedCaloHitList().size());

    if (nOccupiedLayersLhs != nOccupiedLayersRhs)
        return (nOccupiedLayersLhs > nOccupiedLayersRhs);

    const unsigned int layerSpanLhs(pLhs->GetOuterPseudoLayer() - pLhs->GetInnerPseudoLayer());
    const unsigned int layerSpanRhs(pRhs->GetOuterPseudoLayer() - pRhs->GetInnerPseudoLayer());

    if (layerSpanLhs != layerSpanRhs)
        return (layerSpanLhs > layerSpanRhs);

    return (pLhs->GetHadronicEnergy() > pRhs->GetHadronicEnergy());
}


//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::SortByNHits(const Cluster *const pLhs, const Cluster *const pRhs)
{
    const unsigned int nHitsLhs(pLhs->GetNCaloHits());
    const unsigned int nHitsRhs(pRhs->GetNCaloHits());

    if (nHitsLhs != nHitsRhs)
        return (nHitsLhs > nHitsRhs);

    const unsigned int layerSpanLhs(pLhs->GetOuterPseudoLayer() - pLhs->GetInnerPseudoLayer());
    const unsigned int layerSpanRhs(pRhs->GetOuterPseudoLayer() - pRhs->GetInnerPseudoLayer());

    if (layerSpanLhs != layerSpanRhs)
        return (layerSpanLhs > layerSpanRhs);

    return (pLhs->GetHadronicEnergy() > pRhs->GetHadronicEnergy());
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArClusterHelper::TwoDSlidingFitResult::TwoDSlidingFitResult() :
    m_pCluster(NULL),
    m_layerFitHalfWindow(0),
    m_axisIntercept(0.f, 0.f, 0.f),
    m_axisDirection(0.f, 0.f, 0.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArClusterHelper::TwoDSlidingFitResult::LayerFitResult::LayerFitResult(const double l, const double fitT, const double gradient, const double rms) :
    m_l(l),
    m_fitT(fitT),
    m_gradient(gradient),
    m_rms(rms)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

LArClusterHelper::TwoDSlidingFitResult::LayerFitContribution::LayerFitContribution() :
    m_sumT(0.),
    m_sumL(0.),
    m_sumTT(0.),
    m_sumLT(0.),
    m_sumLL(0.),
    m_nPoints(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::TwoDSlidingFitResult::LayerFitContribution::AddPoint(const float l, const float t)
{
    m_sumT += t;
    m_sumL += l;
    m_sumTT += t * t;
    m_sumLT += l * t;
    m_sumLL += l * l;
    ++m_nPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::TwoDSlidingFitResult::GetLocalCoordinates(const CartesianVector &position, float &rL, float &rT) const
{
    const CartesianVector displacement(position - m_axisIntercept);
    const CartesianVector crossProduct(displacement.GetCrossProduct(m_axisDirection));

    rL = displacement.GetDotProduct(m_axisDirection);
    rT = (crossProduct.GetY() < 0.f) ? (-1.f * crossProduct.GetMagnitude()) : crossProduct.GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::TwoDSlidingFitResult::GetGlobalCoordinates(const float rL, const float rT, CartesianVector &position) const
{
    const CartesianVector positiveTDirection(m_axisDirection.GetCrossProduct(CartesianVector(0.f, 1.f, 0.f)));
    position = m_axisIntercept + (m_axisDirection * rL) + (positiveTDirection * rT);
}

//------------------------------------------------------------------------------------------------------------------------------------------

int LArClusterHelper::TwoDSlidingFitResult::GetLayer(const float rL) const
{
    return std::floor(rL / LArPseudoLayerCalculator::GetZPitch());
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::TwoDSlidingFitResult::GetL(const int layer) const
{
    return static_cast<float>(layer) *  LArPseudoLayerCalculator::GetZPitch();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::TwoDSlidingFitResult::GetLocalFitCoordinates(const float x, float &rL, float &rT, int &layer) const
{
    CartesianVector position(0.f, 0.f, 0.f);
    this->GetGlobalFitCoordinates(x, position);
    this->GetLocalCoordinates(position, rL, rT);
    layer = this->GetLayer(rL);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArClusterHelper::TwoDSlidingFitResult::GetGlobalFitCoordinates(const float x, CartesianVector &position) const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    // Find start layer - note this method does not yet work for XZ fits
    if (std::fabs(m_axisDirection.GetX()) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    const float firstL((x - m_axisIntercept.GetX()) / m_axisDirection.GetX());
    const int minLayer(m_layerFitResultMap.begin()->first), maxLayer(m_layerFitResultMap.rbegin()->first);
    const int startLayer(this->GetLayer(firstL));

    if ((startLayer < minLayer) || (startLayer > maxLayer))
        throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

    if (startLayer == minLayer)
    {
        position = this->GetGlobalMinLayerPosition();
        return;
    }

    if (startLayer == maxLayer)
    {
        position = this->GetGlobalMaxLayerPosition();
        return;
    }

    // First layer coordinates
    LayerFitResultMap::const_iterator firstLayerIter(m_layerFitResultMap.end());

    for (int iLayer = startLayer; iLayer < maxLayer; ++iLayer)
    {
        firstLayerIter = m_layerFitResultMap.find(iLayer);

        if (m_layerFitResultMap.end() != firstLayerIter)
            break;
    }

    if (m_layerFitResultMap.end() == firstLayerIter)
    {
std::cout << " first not found " << std::endl;
std::cout << " ->minLayer " << minLayer << " maxLayer " << maxLayer << " startLayer " << startLayer << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }

    const int firstLayer(firstLayerIter->first);
    const float firstLayerL(firstLayerIter->second.GetL());
    const float firstLayerT(firstLayerIter->second.GetFitT());

    CartesianVector firstLayerPosition(0.f, 0.f, 0.f);
    this->GetGlobalCoordinates(firstLayerL, firstLayerT, firstLayerPosition);

    // Sort out directions
    const bool firstIsAheadInX(firstLayerPosition.GetX() > x);
    const bool xIncreasesWithLayers(m_axisDirection.GetX() > 0.f);
    const int increment = ((firstIsAheadInX == xIncreasesWithLayers) ? -1 : +1);

    // Second layer coordinates
    LayerFitResultMap::const_iterator secondLayerIter(m_layerFitResultMap.end());

    for (int iLayer = firstLayer + increment; (iLayer >= minLayer) && (iLayer <= maxLayer); iLayer += increment)
    {
        LayerFitResultMap::const_iterator tempIter = m_layerFitResultMap.find(iLayer);

        if (m_layerFitResultMap.end() == tempIter)
            continue;

        secondLayerIter = tempIter;
        const int secondLayer(secondLayerIter->first);
        const float secondLayerL(secondLayerIter->second.GetL());
        const float secondLayerT(secondLayerIter->second.GetFitT());

        CartesianVector secondLayerPosition(0.f, 0.f, 0.f);
        this->GetGlobalCoordinates(secondLayerL, secondLayerT, secondLayerPosition);

        if ((firstIsAheadInX && (secondLayerPosition.GetX() < x)) || (!firstIsAheadInX && (secondLayerPosition.GetX() > x)))
            break;
    }

    if (m_layerFitResultMap.end() == secondLayerIter)
    {
std::cout << " second not found " << std::endl;
std::cout << " ->minLayer " << minLayer << " maxLayer " << maxLayer << " startLayer " << startLayer << " firstLayer " << firstLayer << " increment=" << increment << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
    }

    const int secondLayer(secondLayerIter->first);
    const float secondLayerL(secondLayerIter->second.GetL());
    const float secondLayerT(secondLayerIter->second.GetFitT());

    CartesianVector secondLayerPosition(0.f, 0.f, 0.f);
    this->GetGlobalCoordinates(secondLayerL, secondLayerT, secondLayerPosition);

    // Linear interpolation
    const float deltaX(x - firstLayerPosition.GetX()), deltaXLayers(secondLayerPosition.GetX() - firstLayerPosition.GetX());

    if (std::fabs(deltaXLayers) < std::numeric_limits<float>::epsilon())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    position = firstLayerPosition + (secondLayerPosition - firstLayerPosition) * (deltaX / deltaXLayers);
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArClusterHelper::TwoDSlidingFitResult::GetGlobalMinLayerPosition() const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    LayerFitResultMap::const_iterator iter = m_layerFitResultMap.begin();
    CartesianVector position(0.f, 0.f, 0.f);
    this->GetGlobalCoordinates(iter->second.GetL(), iter->second.GetFitT(), position);
    return position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector LArClusterHelper::TwoDSlidingFitResult::GetGlobalMaxLayerPosition() const
{
    if (m_layerFitResultMap.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    LayerFitResultMap::const_reverse_iterator iter = m_layerFitResultMap.rbegin();
    CartesianVector position(0.f, 0.f, 0.f);
    this->GetGlobalCoordinates(iter->second.GetL(), iter->second.GetFitT(), position);
    return position;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LArClusterHelper::TwoDSlidingFitResult::IsMultivaluedInX() const
{
  CartesianVector previousPosition(0.f, 0.f, 0.f);
    unsigned int nSteps(0), nPositiveSteps(0), nNegativeSteps(0), nUnchangedSteps(0);

    for (LayerFitResultMap::const_iterator iter = m_layerFitResultMap.begin(), iterEnd = m_layerFitResultMap.end(); iter != iterEnd; ++iter)
    {
        CartesianVector position(0.f, 0.f, 0.f);
        this->GetGlobalCoordinates(iter->second.GetL(), iter->second.GetFitT(), position);

        const CartesianVector delta(position - previousPosition);
        previousPosition = position;
        ++nSteps;

        if (std::fabs(delta.GetX()) < std::fabs(delta.GetZ()) * LArClusterHelper::m_multiValuedTanThetaCut)
        {
            ++nUnchangedSteps;
        }
        else if (delta.GetX() > 0.f)
        {
            ++nPositiveSteps;
        }
        else
        {
            ++nNegativeSteps;
        }
    }

    if (0 == nSteps)
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    std::cout << " nSteps " << nSteps << " nUnchangedSteps " << nUnchangedSteps << " nPositiveSteps " << nPositiveSteps << " nNegativeSteps " << nNegativeSteps << std::endl;    

    if (static_cast<float>(nPositiveSteps) / static_cast<float>(nSteps) > 0.5 || static_cast<float>(nNegativeSteps) / static_cast<float>(nSteps))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArClusterHelper::TwoDSlidingFitResult::GetSlidingFitWidth() const
{
    FloatVector residuals;
    const OrderedCaloHitList &orderedCaloHitList(m_pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            float rL(0.f), rT(0.f);
            this->GetLocalCoordinates((*hitIter)->GetPositionVector(), rL, rT);
            const int layer(this->GetLayer(rL));

            LayerFitResultMap::const_iterator fitResultIter = m_layerFitResultMap.find(layer);

            if (m_layerFitResultMap.end() == fitResultIter)
                continue;

            const double fitT(fitResultIter->second.GetFitT());
            const double gradient(fitResultIter->second.GetGradient());
            const double residualSquared((fitT - rT) * (fitT - rT) / (1. + gradient * gradient)); // angular correction (note: this is cheating!)
            residuals.push_back(residualSquared);
        }
    }

    if (residuals.empty())
        throw StatusCodeException(STATUS_CODE_NOT_INITIALIZED);

    std::sort(residuals.begin(), residuals.end());
    static const float m_trackResidualQuantile(0.8f);
    const float theQuantile(residuals[m_trackResidualQuantile * residuals.size()]);

    return std::sqrt(theQuantile);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArClusterHelper::TwoDSlidingFitResult::FindLargestScatter(unsigned int &largestScatterLayer) const
{
    // Bail out if track is too short
    const unsigned int nFitLayers(m_layerFitResultMap.size());

    if ((nFitLayers <= 2) || (nFitLayers <= 2 * m_layerFitHalfWindow))
        return STATUS_CODE_NOT_FOUND;

    const int minLayer(m_layerFitResultMap.begin()->first), maxLayer(m_layerFitResultMap.rbegin()->first);
    const int nLayersSpanned(1 + maxLayer - minLayer);
    const int layerFitHalfWindow(m_layerFitHalfWindow);

    // Find point of largest scatter
    double splitCosTheta(m_minCosScatteringAngle);
    LayerFitResultMap::const_iterator splitLayerIter(m_layerFitResultMap.end());

    for (LayerFitResultMap::const_iterator iter1 = m_layerFitResultMap.begin(); iter1 != m_layerFitResultMap.end(); ++iter1)
    {
        if (iter1->first - minLayer >= nLayersSpanned - 2 * layerFitHalfWindow)
            break;

        LayerFitResultMap::const_iterator iter2 = m_layerFitResultMap.find(iter1->first + 2 * layerFitHalfWindow);

        if (m_layerFitResultMap.end() == iter2)
            continue;

        const double r1(iter1->second.GetRms()), r2(iter2->second.GetRms());
        const double m1(iter1->second.GetGradient()), m2(iter2->second.GetGradient());
        const double cosTheta = (1. + m1 * m2) / (std::sqrt(1. + m1 * m1) * std::sqrt(1. + m2 * m2));

        if ((r1 < LArClusterHelper::m_trackFitMaxRms) && (r2 < LArClusterHelper::m_trackFitMaxRms) && (cosTheta < splitCosTheta))
        {
            splitCosTheta = cosTheta;

            // Find occupied layer at centre of the kink
            for (int iLayer = iter1->first + layerFitHalfWindow; iLayer < maxLayer; ++iLayer)
            {
                splitLayerIter = m_layerFitResultMap.find(iLayer);

                if (m_layerFitResultMap.end() != splitLayerIter)
                    break;
            }
        }
    }

    if (m_layerFitResultMap.end() == splitLayerIter)
        return STATUS_CODE_NOT_FOUND;

    CartesianVector splitPosition(0.f, 0.f, 0.f);
    this->GetGlobalCoordinates(splitLayerIter->second.GetL(), splitLayerIter->second.GetFitT(), splitPosition);
    largestScatterLayer = GeometryHelper::GetPseudoLayer(splitPosition);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int LArClusterHelper::m_layerFitHalfWindow = 20;
float LArClusterHelper::m_trackFitMaxRms = 0.25f; // cm
float LArClusterHelper::m_minCosScatteringAngle = std::cos(M_PI * 20.f / 180.f); // radians
float LArClusterHelper::m_multiValuedTanThetaCut = 0.1f;
float LArClusterHelper::m_multiValuedStepFractionCut = 0.5f;

StatusCode LArClusterHelper::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LayerFitHalfWindow", m_layerFitHalfWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "TrackFitMaxRms", m_trackFitMaxRms));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCosScatteringAngle", m_minCosScatteringAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MultiValuedTanThetaCut", m_multiValuedTanThetaCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MultiValuedStepFractionCut", m_multiValuedStepFractionCut));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar
