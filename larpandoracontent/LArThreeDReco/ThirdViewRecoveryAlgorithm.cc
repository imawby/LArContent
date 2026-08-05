
/**
 *
 *  @brief  Implementation of the cheating cluster creation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArThreeDReco/ThirdViewRecoveryAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ThirdViewRecoveryAlgorithm::ThirdViewRecoveryAlgorithm() :
    m_minNCaloHits(5),
    m_slidingFitWindow(20),
    m_matchedClusterMaxSep(1.f),
    m_gapTolerance(0.f),
    m_maxMatchedHitSep(1.f),    
    m_minMatchedFrac(0.8f),
    m_maxRecoveryIterations(5),
    m_recoveryMaxTransSep(1.f),
    m_keepMaxTransSep(2.f),
    m_matchedXRange(0.5f),
    m_thresholdOverlapFracForCompatibility(0.8f),
    m_thresholdMatchedFracForCompatibility(0.8f),
    m_stepSize(0.5f),
    m_maxChi2ForMatch(1.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThirdViewRecoveryAlgorithm::Run()
{
    const PfoList *pShowerPfos(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_showerPfoListName, pShowerPfos));

    if (!pShowerPfos || pShowerPfos->empty())
        return STATUS_CODE_NOT_INITIALIZED;

    for (const ParticleFlowObject *const pPfo : *pShowerPfos)
    {
        if (LArPfoHelper::GetNViews(pPfo) != 2)
            continue;
            
        this->RecoverThirdView(pPfo);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThirdViewRecoveryAlgorithm::RecoverThirdView(const Pfo *const pPfo)
{
    std::vector<HitType> hitTypes(this->GetViews(pPfo));
    const Cluster *const pCluster1(this->Get2DCluster(pPfo, hitTypes.at(0)));
    const Cluster *const pCluster2(this->Get2DCluster(pPfo, hitTypes.at(1)));

    if ((pCluster1->GetNCaloHits() < m_minNCaloHits) || (pCluster2->GetNCaloHits() < m_minNCaloHits))
        return;
    
    // Project clusters into third/missing view
    float minX(0.f), maxX(0.f);
    CartesianPointVector projection;
    if (this->GetThirdViewProjection(pCluster1, pCluster2, minX, maxX, projection) != STATUS_CODE_SUCCESS)
        return;

    if (projection.empty())
        return;
    
    // Find closest cluster...
    const Cluster *pMatchedCluster(nullptr);
    this->GetMatchedCluster(projection, hitTypes, pMatchedCluster);

    if (!pMatchedCluster || (pMatchedCluster->GetNCaloHits() < m_minNCaloHits))
        return;

    if (!this->DoesClusterMatchProjections(projection, pMatchedCluster))
        return;
    
    // Look for hits to steal
    CaloHitList collectedHits;
    const Pfo *pMatchedPfo(nullptr);
    
    if (pMatchedCluster->IsAvailable())
    {
        this->RecoverHitsFromAvailable(pMatchedCluster, projection, minX, maxX, collectedHits);
    }
    else
    {
        // Get parent pfo
        this->GetParentPfo(pMatchedCluster, pMatchedPfo);
        int nViews(LArPfoHelper::GetNViews(pMatchedPfo));

        if (nViews == 3)
        {
            this->RecoverHitsFromThreeView(pMatchedPfo, pMatchedCluster, projection, hitTypes, minX, maxX, collectedHits);
        }
        else
        {
            this->RecoverHitsFromTwoView(pMatchedPfo, pMatchedCluster, projection, minX, maxX, collectedHits);
        }
    }

    if (collectedHits.empty())
        return;
    
    // Reassign hits
    this->ReassignHits(pPfo, collectedHits, pMatchedCluster);
}


//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ThirdViewRecoveryAlgorithm::GetThirdViewProjection(const Cluster *const pCluster1, const Cluster *const pCluster2,
    float &xMin, float &xMax, CartesianPointVector &projection)
{
    CaloHitList caloHitList1, caloHitList2;
    LArClusterHelper::GetAllHits(pCluster1, caloHitList1);
    LArClusterHelper::GetAllHits(pCluster2, caloHitList2);    
    HitType hitType1(LArClusterHelper::GetClusterHitType(pCluster1)), hitType2(LArClusterHelper::GetClusterHitType(pCluster2));
    projection.reserve(caloHitList1.size() + caloHitList2.size());
    
    try
    {
        const float slidingFitPitch1(LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType1));
        const TwoDSlidingFitResult slidingFit1(pCluster1, m_slidingFitWindow, slidingFitPitch1);
        const float slidingFitPitch2(LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType2));
        const TwoDSlidingFitResult slidingFit2(pCluster2, m_slidingFitWindow, slidingFitPitch2);

        float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
        pCluster1->GetClusterSpanX(xMin1, xMax1);
        pCluster2->GetClusterSpanX(xMin2, xMax2);
        xMin = std::max(xMin1, xMin2);
        xMax = std::min(xMax1, xMax2);

        // Add projections from first cluster
        this->ProjectHitsToThirdViewWithFits(caloHitList1, slidingFit1, slidingFit2, xMin, xMax, projection);
        // Add projections from second cluster
        this->ProjectHitsToThirdViewWithFits(caloHitList2, slidingFit1, slidingFit2, xMin, xMax, projection);        
    }
    catch (...)
    {
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThirdViewRecoveryAlgorithm::ProjectHitsToThirdViewWithFits(const CaloHitList &caloHitList, const TwoDSlidingFitResult &slidingFit1,
    const TwoDSlidingFitResult &slidingFit2, const float xMin, const float xMax, CartesianPointVector &projection)
{
    const HitType hitType1(LArClusterHelper::GetClusterHitType(slidingFit1.GetCluster()));
    const HitType hitType2(LArClusterHelper::GetClusterHitType(slidingFit2.GetCluster()));
    
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const float thisX(pCaloHit->GetPositionVector().GetX());
            
        if ((thisX < xMin) || (thisX > xMax))
            continue;

        CartesianVector pos1(0.f, 0.f, 0.f);
        if (STATUS_CODE_SUCCESS != slidingFit1.GetGlobalFitPositionAtX(thisX, pos1))
            continue;

        CartesianVector pos2(0.f, 0.f, 0.f);
        if (STATUS_CODE_SUCCESS != slidingFit2.GetGlobalFitPositionAtX(thisX, pos2))
            continue;

        const float z(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType2, pos1.GetZ(), pos2.GetZ()));
        projection.emplace_back(thisX, 0.f, z);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
// We can be smarter here...
void ThirdViewRecoveryAlgorithm::GetMatchedCluster(const CartesianPointVector &projection, const std::vector<HitType> &hitTypes, const Cluster *&pClosestCluster)
{
    const ClusterList *pClusterList3(nullptr);
    const std::string clusterListName3(hitTypes.at(2) == TPC_VIEW_U ? m_clusterListNameU : (hitTypes.at(2) == TPC_VIEW_V ? m_clusterListNameV : m_clusterListNameW));

    if (PandoraContentApi::GetList(*this, clusterListName3, pClusterList3) != STATUS_CODE_SUCCESS) { return; }
    if (!pClusterList3) { return;}

    int bestNMatched(0);    
    for (const Cluster *const pCluster : *pClusterList3)
    {
        int nMatched(0);
        for (const CartesianVector &projPos : projection)
        {
            const float distance(LArClusterHelper::GetClosestDistance(projPos, pCluster));
            
            if (distance < m_matchedClusterMaxSep)
                ++nMatched;
        }

        if (nMatched > bestNMatched)
        {
            bestNMatched = nMatched;
            pClosestCluster = pCluster;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ThirdViewRecoveryAlgorithm::DoesClusterMatchProjections(const CartesianPointVector &projections, const Cluster *const pMatchedCluster)
{
    CaloHitList matchedHitList;
    pMatchedCluster->GetOrderedCaloHitList().FillCaloHitList(matchedHitList);
    HitType matchedHitType(LArClusterHelper::GetClusterHitType(pMatchedCluster));
    
    int nGoodProjPos(0), nInGap(0);
    for (const CartesianVector &projPos : projections)
    {
        if (LArGeometryHelper::IsInGap(this->GetPandora(), projPos, matchedHitType, m_gapTolerance))
        {
            ++nInGap;
            continue;
        }
        
        const float closestDist(LArClusterHelper::GetClosestDistance(projPos, matchedHitList));

        if (closestDist < m_maxMatchedHitSep)
            ++nGoodProjPos;
    }

    const int nValidProj(projections.size() - nInGap);
    if (nValidProj == 0) { return false; }
    const float matchedFrac(float(nGoodProjPos) / float(nValidProj));
    
    return (matchedFrac > m_minMatchedFrac);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThirdViewRecoveryAlgorithm::RecoverHitsFromAvailable(const Cluster *const pMatchedCluster, const CartesianPointVector &projection,
    const float minX, const float maxX, CaloHitList &collectedHits)
{
    this->RecoverHitsWithoutMatchedClusterFit(pMatchedCluster, projection, minX, maxX, collectedHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThirdViewRecoveryAlgorithm::RecoverHitsWithoutMatchedClusterFit(const Cluster *const pMatchedCluster, const CartesianPointVector &projections,
    const float minX, const float maxX, CaloHitList &foundCaloHitList)
{
    foundCaloHitList.clear();
    
    CaloHitList matchedHitList;
    LArClusterHelper::GetAllHits(pMatchedCluster, matchedHitList);

    CartesianPointVector fitPositions;
    for (const CartesianVector &projPos : projections)
        fitPositions.emplace_back(projPos);

    bool found(true);
    while(found)
    {
        found = false;
        
        try
        {
            const float slidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pMatchedCluster)));
            const TwoDSlidingFitResult slidingFitResult(&fitPositions, m_slidingFitWindow, slidingFitPitch);

            for (const CaloHit *const pCaloHit : matchedHitList)
            {
                if ((pCaloHit->GetPositionVector().GetX() < minX) || (pCaloHit->GetPositionVector().GetX() > maxX))
                    continue;
                
                if (std::find(foundCaloHitList.begin(), foundCaloHitList.end(), pCaloHit) != foundCaloHitList.end())
                    continue;
                
                float l(0.f), t(0.f);
                slidingFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), l, t);
                
                if (std::fabs(t) > m_recoveryMaxTransSep)
                    continue;

                found = true;
                foundCaloHitList.emplace_back(pCaloHit);
                fitPositions.emplace_back(pCaloHit->GetPositionVector());
            }
        }
        catch(...)
        {
            break;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void ThirdViewRecoveryAlgorithm::GetParentPfo(const Cluster *const pMatchedCluster, const Pfo *&pMatchedPfo)
{
    const PfoList *pTrackPfos(nullptr), *pShowerPfos(nullptr);
    PandoraContentApi::GetList(*this, m_trackPfoListName, pTrackPfos);
    PandoraContentApi::GetList(*this, m_showerPfoListName, pShowerPfos);
    const HitType hitType3(LArClusterHelper::GetClusterHitType(pMatchedCluster));    

    bool found(false);     
    for (const PfoList *const pPfoList : {pTrackPfos, pShowerPfos})
    {
        if (!pPfoList || pPfoList->empty())
            continue;

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            ClusterList clusters;
            LArPfoHelper::GetClusters(pPfo, hitType3, clusters);

            if (clusters.empty())
                continue;
            
            if (clusters.front() == pMatchedCluster)
            {
                pMatchedPfo = pPfo;
                found = true;
                break;
            }
        }

        if (found) { break; }
    }   
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThirdViewRecoveryAlgorithm::RecoverHitsFromThreeView(const Pfo *const pMatchedPfo, const Cluster *const pMatchedCluster, const CartesianPointVector &projection,
    const std::vector<HitType> &hitTypes, const float minX, const float maxX, CaloHitList &collectedHits)
{
    // Get clusters of the matched pfo in the views of the pfo to recover   
    ClusterList matchedClusters1, matchedClusters2;
    LArPfoHelper::GetClusters(pMatchedPfo, hitTypes.at(0), matchedClusters1);
    LArPfoHelper::GetClusters(pMatchedPfo, hitTypes.at(1), matchedClusters2);
    if (matchedClusters1.empty() || matchedClusters2.empty()) { return; }
    const Cluster *const pMatched1(matchedClusters1.front()), *const pMatched2(matchedClusters2.front());

    // Get matched projections for matched pfo
    float matchedMinX(0.f), matchedMaxX(0.f);
    CartesianPointVector matchedProjection;
    if (this->GetThirdViewProjection(pMatched1, pMatched2, matchedMinX, matchedMaxX, matchedProjection) != STATUS_CODE_SUCCESS)
        return;

    if (matchedProjection.empty())
        return;

    this->RecoverHitsWithMatchedClusterFit(pMatchedCluster, matchedProjection, projection, minX, maxX, collectedHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void ThirdViewRecoveryAlgorithm::RecoverHitsWithMatchedClusterFit(const Cluster *const pMatchedCluster, const CartesianPointVector &matchedProjections,
    const CartesianPointVector &projections, const float minX, const float maxX, CaloHitList &foundCaloHitList)
{
    foundCaloHitList.clear();
    
    CaloHitList matchedHitList;
    LArClusterHelper::GetAllHits(pMatchedCluster, matchedHitList);

    CartesianPointVector fitPositions;
    for (const CartesianVector &projPos : projections)
        fitPositions.emplace_back(projPos);

    // Matched pfo fit in third (recovery) view, using information from the other two views
    const float matchedSlidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pMatchedCluster)));
    const TwoDSlidingFitResult matchedSlidingFitResult(&matchedProjections, m_slidingFitWindow, matchedSlidingFitPitch);
    const float matchedMinL(matchedSlidingFitResult.GetL(matchedSlidingFitResult.GetMinLayer()));
    const float matchedMaxL(matchedSlidingFitResult.GetL(matchedSlidingFitResult.GetMaxLayer()));    

    int iteration(0);
    bool found(true);
    
    while(found && (iteration++ < m_maxRecoveryIterations))
    {
        found = false;
        
        try
        {
            const float slidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pMatchedCluster)));
            const TwoDSlidingFitResult slidingFitResult(&fitPositions, m_slidingFitWindow, slidingFitPitch);

            for (const CaloHit *const pCaloHit : matchedHitList)
            {
                if ((pCaloHit->GetPositionVector().GetX() < minX) || (pCaloHit->GetPositionVector().GetX() > maxX))
                    continue;
                
                if (std::find(foundCaloHitList.begin(), foundCaloHitList.end(), pCaloHit) != foundCaloHitList.end())
                    continue;
                
                float l(0.f), t(0.f);
                slidingFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), l, t);

                float matchedL(0.f), matchedT(0.f);
                matchedSlidingFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), matchedL, matchedT);                
                
                // Does hit fit recovery pfo fit better than matched pfo fit?
                if ((std::fabs(t) > m_recoveryMaxTransSep) && (matchedL > matchedMinL) && (matchedL < matchedMaxL) && (std::fabs(matchedT) < m_keepMaxTransSep))
                    continue;

                found = true;
                foundCaloHitList.emplace_back(pCaloHit);
                fitPositions.emplace_back(pCaloHit->GetPositionVector());
            }
        }
        catch(...)
        {
            break;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThirdViewRecoveryAlgorithm::RecoverHitsFromTwoView(const Pfo *const pMatchedPfo, const Cluster *const pMatchedCluster, const CartesianPointVector &projection,
    const float minX, const float maxX, CaloHitList &collectedHits)
{
    this->RecoverHitsWithoutMatchedClusterFit(pMatchedCluster, projection, minX, maxX, collectedHits);

    // We can likely steal these hits, but we need to rule out that the track & track/shower genuinley overlap in this view
    // We only have two views for each cluster, so we are going to:
    // 1. Identify the 'other' view of the matched two-view pfo and project into third view
    std::vector<HitType> matchedHitTypes(this->GetViews(pMatchedPfo));
    const Cluster *pOtherMatchedCluster(this->Get2DCluster(pMatchedPfo, matchedHitTypes.at(0)) == pMatchedCluster ?
        this->Get2DCluster(pMatchedPfo, matchedHitTypes.at(1)) : this->Get2DCluster(pMatchedPfo, matchedHitTypes.at(0)));
    
    // 2. Project the 'recovered' hits with the other view of the matched pfo
    CartesianPointVector matchedProjections;
    CaloHitList matchedOtherHits;
    LArClusterHelper::GetAllHits(pOtherMatchedCluster, matchedOtherHits);    
    this->ProjectHitsToThirdViewWithHits(matchedOtherHits, collectedHits, matchedProjections);

    // 3. Collect any third view hits
    CaloHitList matchedCollectedHits;
    float matchedMinX(std::numeric_limits<float>::max());
    float matchedMaxX(std::numeric_limits<float>::lowest());
    this->GetMatchedHitsFromView(matchedProjections, matchedHitTypes.at(2), matchedMinX, matchedMaxX, matchedCollectedHits);

    // 4. Investigate collected hits to check validity of match for genuine overlap
    // If we map on to nothing, then great!
    if (matchedCollectedHits.empty())
        return;

    // Bad span? Leave with collected hits!
    const float overlapMin(std::max(matchedMinX, minX));
    const float overlapMax(std::min(matchedMaxX, maxX));
    const float overlap(overlapMax - overlapMin);
    const float overlapFrac(overlap / (maxX - minX));

    if (overlapFrac < m_thresholdOverlapFracForCompatibility) { return; }

    // If good  overlap, test projection...
    const float matchedFraction(this->CalculateThreeViewMatchFraction(pMatchedCluster, pOtherMatchedCluster, matchedHitTypes, matchedCollectedHits, overlapMin, overlapMax));
    
    // Bad match? Leave with collected hits!
    if (matchedFraction < m_thresholdMatchedFracForCompatibility) { return; }
    
    // Seems to be good match.. lets leave the hits in this pfo
    collectedHits.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThirdViewRecoveryAlgorithm::ProjectHitsToThirdViewWithHits(const CaloHitList &smallCaloHitList, const CaloHitList &bigCaloHitList, CartesianPointVector &projections)
{
    std::vector<std::pair<float, const CaloHit*>> bigCaloHitsX;
    bigCaloHitsX.reserve(bigCaloHitList.size());
    
    for (const auto* hit : bigCaloHitList)
        bigCaloHitsX.emplace_back(hit->GetPositionVector().GetX(), hit);
    std::sort(bigCaloHitsX.begin(), bigCaloHitsX.end(), [](auto &a, auto &b) { return a.first < b.first; });

    for (const CaloHit *const pSmallHit : smallCaloHitList)
    {
        // lower bound: first element >= x1 - m_matchedXRange
        auto lower = std::lower_bound(bigCaloHitsX.begin(), bigCaloHitsX.end(), pSmallHit->GetPositionVector().GetX() - m_matchedXRange,
        [](const auto &pair, float value)
        {
            return pair.first < value;
        });

        // scan until we exceed x1 + m_matchedXRange
        for (auto it = lower; it != bigCaloHitsX.end(); ++it)
        {            
            if (it->first > pSmallHit->GetPositionVector().GetX() + m_matchedXRange) break;

            const CaloHit *const pBigHit(it->second);
            const float z(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), pSmallHit->GetHitType(), pBigHit->GetHitType(),
                pSmallHit->GetPositionVector().GetZ(), pBigHit->GetPositionVector().GetZ()));
            
            projections.emplace_back((pBigHit->GetPositionVector().GetX() + pSmallHit->GetPositionVector().GetX()) * 0.5, 0.f, z);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThirdViewRecoveryAlgorithm::GetMatchedHitsFromView(const CartesianPointVector &projections, const HitType targetView,
    float &minX, float &maxX, CaloHitList &collectedHits)
{
    const CaloHitList *pEventViewHits(nullptr);
    std::string viewCaloHitListName(targetView == TPC_VIEW_U ? m_caloHitListNameU : (targetView == TPC_VIEW_V ? m_caloHitListNameV : m_caloHitListNameW));
    PandoraContentApi::GetList(*this, viewCaloHitListName, pEventViewHits);
    if (!pEventViewHits) { return; }
    CaloHitVector eventHitVec(pEventViewHits->begin(), pEventViewHits->end());
    std::sort(eventHitVec.begin(), eventHitVec.end(), LArClusterHelper::SortHitsByPositionInX);
    
    for (const CartesianVector &projection : projections)
    {
        // lower bound: first element >= x - m_matchedXRange
        auto lower = std::lower_bound(eventHitVec.begin(), eventHitVec.end(), projection.GetX() - m_matchedXRange,
        [](const auto &caloHit, float value)
        {
            return caloHit->GetPositionVector().GetX() < value;
        });
        
        // scan until we exceed x + m_matchedXRange
        for (auto it = lower; it != eventHitVec.end(); ++it)
        {
            if ((*it)->GetPositionVector().GetX() > projection.GetX() + m_matchedXRange) break;
            const CaloHit* hit = (*it);

            if (std::fabs(hit->GetPositionVector().GetZ() - projection.GetZ()) < m_matchedXRange)
            {
                collectedHits.emplace_back(hit);                
                minX = std::min(minX, hit->GetPositionVector().GetX());
                maxX = std::max(maxX, hit->GetPositionVector().GetX());
            }
        }
    }
}
    

//------------------------------------------------------------------------------------------------------------------------------------------

float ThirdViewRecoveryAlgorithm::CalculateThreeViewMatchFraction(const Cluster *const pMatchedCluster, const Cluster *const pOtherMatchedCluster,
    const std::vector<HitType> &matchedHitTypes, const CaloHitList &matchedCollectedHits, const float overlapMin, const float overlapMax)
{
    HitType hitType1(LArClusterHelper::GetClusterHitType(pMatchedCluster)), hitType2(LArClusterHelper::GetClusterHitType(pOtherMatchedCluster)), hitType3(matchedHitTypes.at(2));
    const float slidingFitPitch1(LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType1));
    const TwoDSlidingFitResult slidingFitResult1(pMatchedCluster, m_slidingFitWindow, slidingFitPitch1);
    const float slidingFitPitch2(LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType2));
    const TwoDSlidingFitResult slidingFitResult2(pOtherMatchedCluster, m_slidingFitWindow, slidingFitPitch2);
    CartesianPointVector matchedCollectedPositions;
    for (const CaloHit *const pCaloHit : matchedCollectedHits)
        matchedCollectedPositions.push_back(pCaloHit->GetPositionVector());
    const float slidingFitPitch3(LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType3));
    const TwoDSlidingFitResult slidingFitResult3(&matchedCollectedPositions, m_slidingFitWindow, slidingFitPitch3);  
    
    const float overlap(overlapMax - overlapMin);
    int nSamplingPoints(std::floor(overlap / m_stepSize));
    int matchedSamplingPoints(0);

    for (int i=0; i < nSamplingPoints; ++i)
    {
        float thisX(overlapMin + (float(i) * m_stepSize));
        
        CartesianVector pos1(0.f, 0.f, 0.f), pos2(0.f, 0.f, 0.f), pos3(0.f, 0.f, 0.f);
        
        if (STATUS_CODE_SUCCESS != slidingFitResult1.GetGlobalFitPositionAtX(thisX, pos1))
            continue;
            
        if (STATUS_CODE_SUCCESS != slidingFitResult2.GetGlobalFitPositionAtX(thisX, pos2))
            continue;
        
        if (STATUS_CODE_SUCCESS != slidingFitResult3.GetGlobalFitPositionAtX(thisX, pos3))
            continue;        
        
        const float z_12(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType2, pos1.GetZ(), pos2.GetZ()));
        const float dZ_12(fabs(z_12 - pos3.GetZ()));        
        const float z_13(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType3, pos1.GetZ(), pos3.GetZ()));
        const float dZ_13(fabs(z_13 - pos2.GetZ()));        
        const float z_23(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType2, hitType3, pos2.GetZ(), pos3.GetZ()));
        const float dZ_23(fabs(z_23 - pos1.GetZ()));
        const float chi2((dZ_12 + dZ_13 + dZ_23) / 3.f);
        
        
        if (chi2 < m_maxChi2ForMatch)
            ++matchedSamplingPoints;
    }
    
    return (float(matchedSamplingPoints) / float(nSamplingPoints));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThirdViewRecoveryAlgorithm::ReassignHits(const Pfo *const pPfoToRecover, const CaloHitList &collectedHits, const Cluster *const pMatchedCluster)
{
    CaloHitList matchedCaloHitList;
    pMatchedCluster->GetOrderedCaloHitList().FillCaloHitList(matchedCaloHitList);

    CaloHitList caloHitsToRemove(collectedHits), isolatedHitsToRemove;
    this->SplitIntoHitsAndIsolated(pMatchedCluster, caloHitsToRemove, isolatedHitsToRemove);
    
    // If collected all hits, just reassign cluster
    if (matchedCaloHitList.size() == caloHitsToRemove.size())
    {
        if (!pMatchedCluster->IsAvailable())
        {
            // Get parent pfo
            const Pfo *pMatchedPfo(nullptr);
            this->GetParentPfo(pMatchedCluster, pMatchedPfo);
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromPfo(*this, pMatchedPfo, pMatchedCluster));
        }
        
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfoToRecover, pMatchedCluster));
    }
    else
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pMatchedCluster));
        std::string clusterListName(hitType == TPC_VIEW_U ? m_clusterListNameU : (hitType == TPC_VIEW_V ? m_clusterListNameV : m_clusterListNameW));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));

        const Cluster *pNewCluster(nullptr);
        for (const CaloHit *const pCaloHit : collectedHits)
        {
            if (std::find(isolatedHitsToRemove.begin(), isolatedHitsToRemove.end(), pCaloHit) != isolatedHitsToRemove.end())
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveIsolatedFromCluster(*this, pMatchedCluster, pCaloHit));
            }
            else
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, pMatchedCluster, pCaloHit));
            }

            if (!pNewCluster)
            {
                const ClusterList *pTemporaryList(nullptr);
                std::string temporaryListName, currentListName;
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, currentListName));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
                    PandoraContentApi::CreateTemporaryListAndSetCurrent<ClusterList>(*this, pTemporaryList, temporaryListName));

                PandoraContentApi::Cluster::Parameters parameters;
                parameters.m_caloHitList.push_back(pCaloHit);

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewCluster));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, temporaryListName, currentListName));
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, currentListName));
            }
            else
            {
                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, pNewCluster, pCaloHit));
            }
        }

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToPfo(*this, pPfoToRecover, pNewCluster));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void ThirdViewRecoveryAlgorithm::SplitIntoHitsAndIsolated(const Cluster *const pMatchedCluster, CaloHitList &collectedHits, CaloHitList &isolatedCollectedHits)
{
    const CaloHitList &isolatedHits(pMatchedCluster->GetIsolatedCaloHitList());
    
    CaloHitList temp(collectedHits);
    collectedHits.clear();

    for (const CaloHit *const pCaloHit : temp)
    {
        if (std::find(isolatedHits.begin(), isolatedHits.end(), pCaloHit) != isolatedHits.end())
            isolatedCollectedHits.emplace_back(pCaloHit);
        else
            collectedHits.emplace_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster* ThirdViewRecoveryAlgorithm::Get2DCluster(const Pfo *const pPfo, const HitType &hitType)
{
    ClusterList clusters;
    LArPfoHelper::GetClusters(pPfo, hitType, clusters);

    if (clusters.size() == 0)
        throw StatusCodeException(STATUS_CODE_FAILURE);
    
    return clusters.front();
}


//------------------------------------------------------------------------------------------------------------------------------------------    

std::vector<HitType> ThirdViewRecoveryAlgorithm::GetViews(const Pfo *const pPfo)
{
    HitType hitType1(TPC_3D), hitType2(TPC_3D), hitType3(TPC_3D); // is there something better to use?
    
    for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        HitType &hitTypeToSet(hitType1 == TPC_3D ? hitType1 : hitType2);

        ClusterList clusters;        
        LArPfoHelper::GetClusters(pPfo, hitType, clusters);
        int nClusters(clusters.size());
        
        if (nClusters == 0)
            hitType3 = hitType;
        else if (nClusters == 1)
            hitTypeToSet = hitType;
        else
            throw StatusCodeException(STATUS_CODE_FAILURE);
    }
    
    return std::vector<HitType>({hitType1, hitType2, hitType3});
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------    

StatusCode ThirdViewRecoveryAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "TrackPfoListName", m_trackPfoListName));
    if (m_trackPfoListName.empty()) { m_trackPfoListName = "TrackParticles3D"; }
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
    if (m_showerPfoListName.empty()) { m_showerPfoListName = "ShowerParticles3D"; }
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameU", m_caloHitListNameU));
    if (m_caloHitListNameU.empty()) { m_caloHitListNameU = "CaloHitListU"; }
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameV", m_caloHitListNameV));
    if (m_caloHitListNameV.empty()) { m_caloHitListNameV = "CaloHitListV"; }
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameW", m_caloHitListNameW));
    if (m_caloHitListNameW.empty()) { m_caloHitListNameW = "CaloHitListW"; }    
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameU", m_clusterListNameU));
    if (m_clusterListNameU.empty()) { m_clusterListNameU = "ClustersU"; }
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameV", m_clusterListNameV));
    if (m_clusterListNameV.empty()) { m_clusterListNameV = "ClustersV"; }
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameW", m_clusterListNameW));
    if (m_clusterListNameW.empty()) { m_clusterListNameW = "ClustersW"; }


    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinNCaloHits", m_minNCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SlidingFitWindow", m_slidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MatchedClusterMaxSep", m_matchedClusterMaxSep));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "GapTolerance", m_gapTolerance));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxMatchedHitSep", m_maxMatchedHitSep));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxRecoveryIterations", m_maxRecoveryIterations));
    
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "RecoveryMaxTransSep", m_recoveryMaxTransSep));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "KeepMaxTransSep", m_keepMaxTransSep));    

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MatchedXRange", m_matchedXRange));


    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ThresholdOverlapFracForCompatibility", m_thresholdOverlapFracForCompatibility));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ThresholdMatchedFracForCompatibility", m_thresholdMatchedFracForCompatibility));

    
    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
