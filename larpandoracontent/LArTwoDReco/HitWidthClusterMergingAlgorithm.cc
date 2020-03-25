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
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

using namespace pandora;

namespace lar_content
{

HitWidthClusterMergingAlgorithm::HitWidthClusterMergingAlgorithm() :
  m_clusterListName(),
  m_maxConstituentHitWidth(0.5),
  m_hitWidthScalingFactor(1.0),
  m_useSlidingLinearFit(false),
  m_layerFitHalfWindow(20),
  m_fittingWeight(10),
  m_minClusterWeight(0.5),      
  m_maxXMergeDistance(5),     
  m_maxZMergeDistance(2),     
  m_minMergeCosOpeningAngle(0.97),
  m_minDirectionDeviationCosAngle(0.90)
{
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetListOfCleanClusters(const ClusterList *const pClusterList, ClusterVector &clusterVector) const
{
    LArHitWidthHelper::ClusterToParametersMapStore* pClusterToParametersMapStore = LArHitWidthHelper::ClusterToParametersMapStore::Instance();
    LArHitWidthHelper::ClusterToParametersMap clusterToParametersMap = pClusterToParametersMapStore->GetMap();

    // clear map if already full i.e. from other view clustering
    if(!clusterToParametersMap.empty())
        clusterToParametersMap.clear();

    for (const Cluster *const pCluster : *pClusterList)
    {
        // the original cluster weight, with no hit scaling or hit padding
        if (LArHitWidthHelper::GetOriginalTotalClusterWeight(pCluster) < m_minClusterWeight)
            continue;

        clusterToParametersMap.insert(std::pair<const Cluster*, LArHitWidthHelper::ClusterParameters>(pCluster, LArHitWidthHelper::ClusterParameters(pCluster, m_maxConstituentHitWidth, m_useSlidingLinearFit, m_hitWidthScalingFactor)));
        clusterVector.push_back(pCluster);
    }

    //ORDER BY MAX EXTREMAL X COORDINATE
    std::sort(clusterVector.begin(), clusterVector.end(), LArHitWidthHelper::SortByHigherXExtrema);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::PopulateClusterAssociationMap(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{
    // ATTN this method assumes that clusters have been sorted by extremal x position (low higherXExtrema -> high higherXExtrema)
    for (ClusterVector::const_iterator iterCurrentCluster = clusterVector.begin(); iterCurrentCluster != clusterVector.end(); ++iterCurrentCluster)
    {
        const Cluster *const pCurrentCluster = *iterCurrentCluster;
        const LArHitWidthHelper::ClusterParameters currentClusterParameters(LArHitWidthHelper::GetClusterParameters(pCurrentCluster));

        for (ClusterVector::const_iterator iterTestCluster = iterCurrentCluster; iterTestCluster != clusterVector.end(); ++iterTestCluster)
        {
            if (iterCurrentCluster == iterTestCluster)
                continue;

            const Cluster *const pTestCluster = *iterTestCluster;
            const LArHitWidthHelper::ClusterParameters testClusterParameters(LArHitWidthHelper::GetClusterParameters(pTestCluster));

            if (!this->AreClustersAssociated(currentClusterParameters, testClusterParameters))
                continue;

            clusterAssociationMap[pCurrentCluster].m_forwardAssociations.insert(pTestCluster);
            clusterAssociationMap[pTestCluster].m_backwardAssociations.insert(pCurrentCluster);
        }
    }

    // remove all 'shortcut' routes
    this->CleanupClusterAssociations(clusterVector, clusterAssociationMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::CleanupClusterAssociations(const ClusterVector &clusterVector, ClusterAssociationMap &clusterAssociationMap) const
{
    // Create temporary map so can delete elements whilst still iterating over them
    ClusterAssociationMap tempMap(clusterAssociationMap);

    for (const Cluster *const pCluster : clusterVector)
    {
        const ClusterAssociationMap::const_iterator primaryMapIter = clusterAssociationMap.find(pCluster);

        if (primaryMapIter == clusterAssociationMap.end())
            continue;

        const ClusterSet &primaryForwardAssociations(primaryMapIter->second.m_forwardAssociations);

        // loop through the primary associations of each cluster
        // remove clusters that are present in secondary associations of other primary associated clusters
        for (const Cluster *const pConsideredCluster : primaryForwardAssociations)
        {
            for (const Cluster *const pPrimaryCluster : primaryForwardAssociations)
            {
                if (pConsideredCluster == pPrimaryCluster)
                    continue;

                const ClusterAssociationMap::const_iterator secondaryMapIter = clusterAssociationMap.find(pPrimaryCluster);

                // if primary cluster has no associations (this shouldn't ever be the case)
                if (secondaryMapIter == clusterAssociationMap.end())
                    continue;

                const ClusterSet &secondaryForwardAssociations(secondaryMapIter->second.m_forwardAssociations);

                // if considered cluster is present in the forward associations of any other primary associated cluster remove from primary associations
                if (secondaryForwardAssociations.find(pConsideredCluster) != secondaryForwardAssociations.end())
                {
                    ClusterSet &tempPrimaryForwardAssociations(tempMap.find(pCluster)->second.m_forwardAssociations);
                    const ClusterSet::const_iterator forwardAssociationToRemove(tempPrimaryForwardAssociations.find(pConsideredCluster));

                    // if association has already been removed
                    if(forwardAssociationToRemove == tempPrimaryForwardAssociations.end())
                        continue;

                    ClusterSet &tempPrimaryBackwardAssociations(tempMap.find(pConsideredCluster)->second.m_backwardAssociations);
                    const ClusterSet::const_iterator backwardAssociationToRemove(tempPrimaryBackwardAssociations.find(pCluster));

                    // if association has already been removed
                    if (backwardAssociationToRemove == tempPrimaryBackwardAssociations.end())
                        continue;

                    tempPrimaryForwardAssociations.erase(forwardAssociationToRemove);
                    tempPrimaryBackwardAssociations.erase(backwardAssociationToRemove);
                }
            }
        }
    }
  clusterAssociationMap = tempMap;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HitWidthClusterMergingAlgorithm::AreClustersAssociated(const LArHitWidthHelper::ClusterParameters &currentFitParameters, const LArHitWidthHelper::ClusterParameters &testFitParameters) const
{
    CartesianVector currentClusterDirection(0,0,0);
    CartesianVector testClusterDirection(0,0,0);

    if (m_useSlidingLinearFit)
    {
        const CartesianPointVector currentConstituentHitPositionVector(LArHitWidthHelper::GetConstituentHitPositionVector(currentFitParameters.GetConstituentHitVector()));
        const CartesianPointVector testConstituentHitPositionVector(LArHitWidthHelper::GetConstituentHitPositionVector(testFitParameters.GetConstituentHitVector()));
        const float layerPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        TwoDSlidingFitResult currentFitResult(&currentConstituentHitPositionVector, m_layerFitHalfWindow, layerPitch);
        TwoDSlidingFitResult testFitResult(&testConstituentHitPositionVector, m_layerFitHalfWindow, layerPitch);

        StatusCode currentStatus(currentFitResult.GetGlobalFitDirectionAtX(currentFitParameters.GetHigherXExtrema().GetX(), currentClusterDirection));
        StatusCode testStatus(testFitResult.GetGlobalFitDirectionAtX(testFitParameters.GetLowerXExtrema().GetX(), testClusterDirection));
      
        if(currentStatus != STATUS_CODE_SUCCESS || testStatus != STATUS_CODE_SUCCESS)
            return false;
    }
    else
    {
        currentClusterDirection = GetClusterDirection(currentFitParameters, currentFitParameters.GetHigherXExtrema());
        testClusterDirection = GetClusterDirection(testFitParameters, testFitParameters.GetLowerXExtrema());
    }

    if (testFitParameters.GetLowerXExtrema().GetX() > (currentFitParameters.GetHigherXExtrema().GetX() + m_maxXMergeDistance))
      return false;

    if (testFitParameters.GetLowerXExtrema().GetZ() > (currentFitParameters.GetHigherXExtrema().GetZ() + m_maxZMergeDistance) || testFitParameters.GetLowerXExtrema().GetZ() < (currentFitParameters.GetHigherXExtrema().GetZ() - m_maxZMergeDistance))
      return false;
    
    if (currentClusterDirection.GetCosOpeningAngle(testClusterDirection) < m_minMergeCosOpeningAngle)
      return false;
    
    // check that the new direction is consistent with the old clusters
    if (m_useSlidingLinearFit)
    {
        const CartesianPointVector currentConstituentHitPositionVector(LArHitWidthHelper::GetConstituentHitPositionVector(currentFitParameters.GetConstituentHitVector()));
        const CartesianPointVector testConstituentHitPositionVector(LArHitWidthHelper::GetConstituentHitPositionVector(testFitParameters.GetConstituentHitVector()));
        CartesianPointVector newConstituentHitPositionVector(currentConstituentHitPositionVector);
        
        newConstituentHitPositionVector.insert(newConstituentHitPositionVector.end(), testConstituentHitPositionVector.begin(), testConstituentHitPositionVector.end());
      
        const float layerPitch(LArGeometryHelper::GetWireZPitch(this->GetPandora()));
        TwoDSlidingFitResult newFitResult(&newConstituentHitPositionVector, m_layerFitHalfWindow, layerPitch);
        CartesianVector newClusterDirectionCurrent(0,0,0);
        CartesianVector newClusterDirectionTest(0,0,0);
        StatusCode currentStatus(newFitResult.GetGlobalFitDirectionAtX(currentFitParameters.GetHigherXExtrema().GetX(), newClusterDirectionCurrent));
        StatusCode testStatus(newFitResult.GetGlobalFitDirectionAtX(testFitParameters.GetLowerXExtrema().GetX(), newClusterDirectionTest));
      
        if (currentStatus != STATUS_CODE_SUCCESS || testStatus != STATUS_CODE_SUCCESS)
            return false;

        if (newClusterDirectionCurrent.GetCosOpeningAngle(currentClusterDirection) < m_minDirectionDeviationCosAngle || newClusterDirectionTest.GetCosOpeningAngle(testClusterDirection) < m_minDirectionDeviationCosAngle)
            return false;
    }
    else
    {
        const LArHitWidthHelper::ConstituentHitVector currentConstituentHitVector(currentFitParameters.GetConstituentHitVector());
        const LArHitWidthHelper::ConstituentHitVector testConstituentHitVector(testFitParameters.GetConstituentHitVector());
        LArHitWidthHelper::ConstituentHitVector newConstituentHitVector(currentConstituentHitVector);
        
        newConstituentHitVector.insert(newConstituentHitVector.end(), testConstituentHitVector.begin(), testConstituentHitVector.end());
      
        LArHitWidthHelper::ClusterParameters newParameters(nullptr, currentFitParameters.GetNumCaloHits() + testFitParameters.GetNumCaloHits(), currentFitParameters.GetTotalWeight() + testFitParameters.GetTotalWeight(), newConstituentHitVector, CartesianVector(0,0,0), CartesianVector(0,0,0));
        const CartesianVector midpoint((currentFitParameters.GetHigherXExtrema() + testFitParameters.GetLowerXExtrema())*0.5);
        CartesianVector newClusterDirection(GetClusterDirection(newParameters, midpoint));
        
        if (newClusterDirection.GetCosOpeningAngle(currentClusterDirection) < m_minDirectionDeviationCosAngle || newClusterDirection.GetCosOpeningAngle(testClusterDirection) < m_minDirectionDeviationCosAngle)
            return false;
    }
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HitWidthClusterMergingAlgorithm::IsExtremalCluster(const bool isForward, const Cluster *const pCurrentCluster,  const Cluster *const pTestCluster) const
{
    const LArHitWidthHelper::ClusterParameters currentClusterParameters(LArHitWidthHelper::GetClusterParameters(pCurrentCluster));
    const LArHitWidthHelper::ClusterParameters testClusterParameters(LArHitWidthHelper::GetClusterParameters(pTestCluster));
    CartesianVector currentHigherXExtrema(currentClusterParameters.GetHigherXExtrema()), testHigherXExtrema(currentClusterParameters.GetHigherXExtrema());
    float currentMaxX(currentHigherXExtrema.GetX()), testMaxX(testHigherXExtrema.GetX());

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

CartesianVector HitWidthClusterMergingAlgorithm::GetClusterDirection(const LArHitWidthHelper::ClusterParameters &clusterFitParameters, const CartesianVector &fitReferencePoint) const 
{
    // if cluster composed of one hit, return a transverse line as fit
    if (clusterFitParameters.GetNumCaloHits() == 1)
        return CartesianVector(1, 0, 0);

    // minimise the longitudinal distance in fit (works best for transverse fits)
    CartesianVector lsTransverseClusterFitDirection(0,0,0);
    CartesianVector lsTransverseIntercept(0,0,0);
    float lsTransverseChiSquared(0);

    // minimise the transverse coordinate in fit (works best for longitudinal fits)
    CartesianVector lsLongitudinalClusterFitDirection(0,0,0);
    CartesianVector lsLongitudinalIntercept(0,0,0);
    float lsLongitudinalChiSquared(0);

    GetWeightedGradient(clusterFitParameters, true, lsTransverseClusterFitDirection, lsTransverseIntercept, lsTransverseChiSquared, fitReferencePoint);
    GetWeightedGradient(clusterFitParameters, false, lsLongitudinalClusterFitDirection, lsLongitudinalIntercept, lsLongitudinalChiSquared, fitReferencePoint);
    
    // return fit with the lowest chi-squared
    if (lsTransverseChiSquared < lsLongitudinalChiSquared)
        return lsTransverseClusterFitDirection;

    return lsLongitudinalClusterFitDirection;
}
  
//------------------------------------------------------------------------------------------------------------------------------------------

void HitWidthClusterMergingAlgorithm::GetWeightedGradient(const LArHitWidthHelper::ClusterParameters &clusterFitParameters, bool isTransverse, CartesianVector &direction, CartesianVector &zIntercept, float &chiSquared, const CartesianVector &fitReferencePoint) const
{
    // cannot make a longitudinal fit to a single hit cluster
    if (!isTransverse && clusterFitParameters.GetNumCaloHits() == 1) 
    {
        std::cout << "WARNING - CANNOT MAKE LONGITUDINAL FIT TO SINGLE HIT CLUSTER" << std::endl;
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
    }

    float weightSum(0);
    float weightedXSum(0);
    float weightedZSum(0);
    bool isXConstant(true);
    bool isZConstant(true);

    LArHitWidthHelper::ConstituentHitVector constituentHitVector(clusterFitParameters.GetConstituentHitVector());

    // sort hits with respect to their distance to the fitReferencePoint (closest -> furthest)
    std::sort(constituentHitVector.begin(), constituentHitVector.end(), LArHitWidthHelper::ConstituentHit::SortByDistanceToPoint(fitReferencePoint));

    // calculate weightedXMean and weightedZMean for fit
    for (const LArHitWidthHelper::ConstituentHit &constituentHit : constituentHitVector) 
    {
        const CartesianVector hitPosition = constituentHit.GetPositionVector();
        const float hitWeight = constituentHit.GetHitWidth();

        if (std::fabs(clusterFitParameters.GetConstituentHitVector().begin()->GetPositionVector().GetX() - hitPosition.GetX()) > std::numeric_limits<float>::epsilon())
            isXConstant = false;

        if (std::fabs(clusterFitParameters.GetConstituentHitVector().begin()->GetPositionVector().GetZ() - hitPosition.GetZ()) > std::numeric_limits<float>::epsilon())
            isZConstant = false;

        weightedXSum += hitPosition.GetX() * hitWeight;
        weightedZSum += hitPosition.GetZ() * hitWeight;
        weightSum += hitWeight;

        // do not exceed the specified hit weight in the fit
        if(weightSum > m_fittingWeight)
            break;
    }

    // return vertical fit for a cluster with constant x
    if (isXConstant) 
    {
        direction = CartesianVector(0, 0, 1);
        zIntercept = CartesianVector(0, 0, 0);
        chiSquared = 0;
        return;
    }

    // return horizontal fit for a cluster with constant z
    if (isZConstant) 
    {
        direction = CartesianVector(1, 0, 0);
        zIntercept = CartesianVector(0, 0, clusterFitParameters.GetConstituentHitVector().begin()->GetPositionVector().GetZ());
        chiSquared = 0;
        return;
    }

    float weightedXMean(weightedXSum/weightSum);
    float weightedZMean(weightedZSum/weightSum);
    
    float numerator(0);
    float denominator(0);

    float weightCount(0);
    for (const LArHitWidthHelper::ConstituentHit &constituentHit : constituentHitVector) 
    {
        const CartesianVector hitPosition = constituentHit.GetPositionVector();
        const float hitWeight = constituentHit.GetHitWidth();
        weightCount += hitWeight;
        
        numerator += hitWeight * (hitPosition.GetX() - weightedXMean) * (hitPosition.GetZ() - weightedZMean);
        denominator += isTransverse ? hitWeight * pow(hitPosition.GetX() - weightedXMean, 2) : hitWeight * pow(hitPosition.GetZ() - weightedZMean, 2);

        // do not exceed the specified hit weight in the fit
        if (weightCount > m_fittingWeight)
            break;
    }

    float gradient = numerator/denominator;
    
    float intercept(0);
    intercept = isTransverse ? weightedZMean - gradient * weightedXMean : weightedXMean - gradient * weightedZMean;

    float chi(0);
    weightCount = 0;
    for (const LArHitWidthHelper::ConstituentHit &constituentHit : constituentHitVector) 
    {
        const CartesianVector hitPosition = constituentHit.GetPositionVector();
        const float hitWeight = constituentHit.GetHitWidth();
        weightCount += hitWeight;
        
        chi += isTransverse ? hitWeight * pow(hitPosition.GetZ() - intercept - gradient * hitPosition.GetX(), 2) : hitWeight * pow(hitPosition.GetX() - intercept - gradient * hitPosition.GetZ(), 2);

        // do not exceed the specified hit weight in the fit
        if (weightCount > m_fittingWeight)
            break;
    }

    // change coordinates to z=mx+c fit and normalise
    direction = isTransverse ? CartesianVector(1.0, 0, gradient).GetUnitVector() : CartesianVector(gradient, 0, 1.0).GetUnitVector();
    zIntercept = isTransverse ? CartesianVector(0, 0, intercept) : CartesianVector(0, 0, -intercept/gradient);
    chiSquared = chi;

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitWidthClusterMergingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "UseSlidingLinearFit", m_useSlidingLinearFit));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "FittingWeight", m_fittingWeight));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "LayerFitHalfWindow", m_layerFitHalfWindow));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinClusterWeight", m_minClusterWeight));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxXMergeDistance", m_maxXMergeDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxZMergeDistance", m_maxZMergeDistance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinMergeCosOpeningAngle", m_minMergeCosOpeningAngle));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinDirectionDeviationCosAngle", m_minDirectionDeviationCosAngle));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxConstituentHitWidth", m_maxConstituentHitWidth));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, 
        "HitWidthScalingFactor", m_hitWidthScalingFactor));

    return ClusterAssociationAlgorithm::ReadSettings(xmlHandle);
}
  
} // namespace lar_content