/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the two dimensional sliding fit splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"
#include "larpandoracontent/LArHelpers/LArPcaHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArUtility/KalmanFilter.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

KalmanSplittingAlgorithm::KalmanFit::KalmanFit(const KalmanFilter2D &kalmanFilter, const int currentWireID) :
    m_kalmanFilter(kalmanFilter),
    m_positions(CartesianPointVector()),
    m_directions(CartesianPointVector()),
    m_currentWireID(currentWireID)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanSplittingAlgorithm::KalmanFit::SaveStep(const int currentWireID)
{
    const CartesianVector kalmanPos(m_kalmanFilter.GetPosition()(0), 0.f, m_kalmanFilter.GetPosition()(1));
    const CartesianVector kalmanDir(m_kalmanFilter.GetDirection()(0), 0.f, m_kalmanFilter.GetDirection()(1));

    m_positions.push_back(kalmanPos);
    m_directions.push_back(kalmanDir);
    m_currentWireID = currentWireID;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanSplittingAlgorithm::KalmanFit::AddPositionAndUpdate(const CartesianVector &position)
{
    Eigen::VectorXd eigenXd(2);
    eigenXd << position.GetX(), position.GetZ();
    m_kalmanFilter.Update(eigenXd);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

KalmanSplittingAlgorithm::KalmanSplittingAlgorithm() :
    m_minClusterLength(10.f),
    m_kalmanDelta(1.f),
    m_kalmanProcessVarCoeff(1.f),
    m_kalmanMeasurementVarCoeff(1.f),
    m_minTransSeparation(5.f),
    m_segmentWindows(6),
    m_minDeviation(10.f),
    m_maxSpread(0.1f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanSplittingAlgorithm::DivideCaloHits(const Cluster *const pCluster, CaloHitList &firstHitList, CaloHitList &secondHitList) const
{
    // Do we want to work on this cluster?
    if (!this->IsTargetCluster(pCluster))
        return STATUS_CODE_NOT_FOUND;

    //////////////////////////////////
    // PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));
    //////////////////////////////////

    // Get KalmanFit
    try
    {
        //////////////////////////////////
        // ClusterList visualiseClusters({pCluster});
        // PANDORA_MONITORING_API(VisualizeClusters(this->GetPandora(), &visualiseClusters, "Cluster", BLACK));
        // PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
        //////////////////////////////////

        KalmanFit kalmanFit(this->PerformKalmanFit(pCluster));

        if (kalmanFit.m_positions.empty())
            return STATUS_CODE_NOT_FOUND;

        // Do sliding fit
        const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
        const HitType view(LArClusterHelper::GetClusterHitType(pCluster));
        const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
        const TwoDSlidingFitResult twoDSlidingFitResult(pCluster, 20, pitch);

        // Search for split position
        CartesianVector splitPosition(0.f, 0.f, 0.f);

        if (STATUS_CODE_SUCCESS == this->FindBestSplitPosition(pCluster, kalmanFit, splitPosition))
        {
            return this->DivideCaloHits(twoDSlidingFitResult, splitPosition, firstHitList, secondHitList);
        }
    }
    catch (StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_FAILURE == statusCodeException.GetStatusCode())
            throw statusCodeException;
    }

    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool KalmanSplittingAlgorithm::IsTargetCluster(const Cluster *const pCluster) const
{
    if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
        return false;

    const float m_maxConstituentHitWidth(0.5f), m_hitWidthScalingFactor(1.f), m_minClusterSparseness(0.5f);
    const unsigned int numberOfProposedConstituentHits(LArHitWidthHelper::GetNProposedConstituentHits(
        pCluster, m_maxConstituentHitWidth, m_hitWidthScalingFactor));

    if (numberOfProposedConstituentHits == 0)
        return false;

    // Avoid wide hit width clusters
    // clusterSparseness [0 -> 1] where a higher value indicates sparseness
    const float clusterSparseness(1.f - (static_cast<float>(pCluster->GetNCaloHits()) / static_cast<float>(numberOfProposedConstituentHits)));

    if (clusterSparseness > m_minClusterSparseness)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
KalmanSplittingAlgorithm::KalmanFit KalmanSplittingAlgorithm::PerformKalmanFit(const Cluster *const pCluster) const
{
    // Get CaloHits
    CaloHitList caloHitList;
    LArClusterHelper::GetAllHits(pCluster, caloHitList);

    // Create wireID map to assess ambiguities
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const HitType view(caloHitList.front()->GetHitType());
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
   
    std::map<int, CaloHitList> caloHitWireMap;
    for (const CaloHit *const pCaloHit : caloHitList)
        caloHitWireMap[std::floor(pCaloHit->GetPositionVector().GetZ() / pitch)].push_back(pCaloHit);

    // Initialise Kalman fit
    Eigen::VectorXd init(2);
    int startWireID(-1);

    if (this->GetFitSeed(caloHitWireMap, init, startWireID) != STATUS_CODE_SUCCESS)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const float processVariance{m_kalmanProcessVarCoeff * pitch * pitch};
    const float measurementVariance{m_kalmanMeasurementVarCoeff * pitch * pitch};
    KalmanFit kalmanFit(KalmanFilter2D(m_kalmanDelta, processVariance, measurementVariance, init), startWireID);

    // Now populate fit..
    for (auto &entry : caloHitWireMap)
    {
        if (entry.first <= kalmanFit.m_currentWireID)
            continue;
         
        // Make next kalman step
        kalmanFit.m_kalmanFilter.Predict();
       
        // Find matched cluster position
        CartesianPointVector matchedPositions;
        if (this->FindMatchedClusterPosition(entry.second, kalmanFit, matchedPositions) != STATUS_CODE_SUCCESS)
            continue;

        if (matchedPositions.size() == 1) // If unambiguous, just add in hit
        {
            kalmanFit.SaveStep(entry.first);
            kalmanFit.AddPositionAndUpdate(matchedPositions.front());

            ///////////////////////////////////
            // const CartesianVector bestPos(matchedPositions.front());
            // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &bestPos, "kalmanPos", GREEN, 2));                        
            ///////////////////////////////////
        }
        else if (matchedPositions.size() > 1)  // If ambiguous, follow each route and pick the straightest
        {
            ///////////////////////////////////
            // PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            // for (const CartesianVector &seed : matchedPositions)
            // {
            //     std::cout << "Seed: " << seed << std::endl;
            //     PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &seed, "kalmanPos", RED, 2));
            // }
            // PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            ///////////////////////////////////

            KalmanFitVector pathways;

            for (const CartesianVector &seed : matchedPositions)
            {
                pathways.push_back(kalmanFit);
                pathways.back().SaveStep(entry.first);
                pathways.back().AddPositionAndUpdate(seed);

                // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &seed, "kalmanPos", RED, 2));

                this->FollowRoute(caloHitWireMap, pathways.back());

                // PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            }

            // PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));

            // Find the best path
            int highestNHits(0);
            float straightest(std::numeric_limits<float>::max());

            for (const KalmanFit &pathwayFit : pathways)
            {
                const int nHits(pathwayFit.m_positions.size());
                const float sigma(this->GetSTD(pathwayFit.m_directions));

                // Pick longest path
                bool best(nHits > highestNHits);

                // If lengths are similar, pick straightest path
                if (highestNHits != 0)
                {
                    const float hitRatio(std::fabs(static_cast<float>(nHits - highestNHits)) / highestNHits);

                    if (hitRatio < 0.1)
                        best = (sigma < straightest);
                }

                if (best)
                {
                    highestNHits = nHits;
                    straightest = sigma;
                    kalmanFit = pathwayFit;
                }
            }

            //////////////////////////////////
            // for (const CartesianVector &jam : kalmanFit.m_positions)
            // {
            //     CartesianVector frog(jam);
            //     PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &frog, "CHOSEN", GREEN, 2));
            // }
            // PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
            //////////////////////////////////
        }
    }

    //////////////////////////////////
    // PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    //////////////////////////////////
    
    return kalmanFit;
}
    
//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanSplittingAlgorithm::GetFitSeed(const std::map<int, CaloHitList> &caloHitWireMap, Eigen::VectorXd &init, int &startWireID) const
{
    for (auto &entry : caloHitWireMap)
    {
        if (entry.second.size() == 1)
        {
            const CartesianVector &seedPos(entry.second.front()->GetPositionVector());    
            init << seedPos.GetX(), seedPos.GetZ();
            startWireID = entry.first;
            break;
        }
    }

    if (startWireID < 0)
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanSplittingAlgorithm::FindMatchedClusterPosition(const CaloHitList &caloHitList, const KalmanFit &kalmanFit, 
    CartesianPointVector &matchedPositions) const
{
    const CartesianVector kalmanPos(kalmanFit.m_kalmanFilter.GetPosition()(0), 0.f, kalmanFit.m_kalmanFilter.GetPosition()(1));
    const CartesianVector kalmanDir(kalmanFit.m_kalmanFilter.GetDirection()(0), 0.f, kalmanFit.m_kalmanFilter.GetDirection()(1));

    std::vector<std::pair<CartesianVector, float>> matchedMap;

    for (const CaloHit * const pCaloHit : caloHitList)
    {
        const CartesianVector hitPosition((pCaloHit->GetCellSize1() > 0.5f) ?
            LArHitWidthHelper::GetClosestPointToLine2D(kalmanPos, kalmanDir, pCaloHit) : pCaloHit->GetPositionVector()); 

        const float t(kalmanDir.GetCrossProduct(kalmanPos - hitPosition).GetMagnitude());

        if (t < m_minTransSeparation)
            matchedMap.push_back(std::make_pair(pCaloHit->GetPositionVector(), t));
    }

    std::sort(matchedMap.begin(), matchedMap.end(), [](const std::pair<CartesianVector, float> &lhs, const std::pair<CartesianVector, float> &rhs) 
        {
            return lhs.second < rhs.second;
        }
    );

    for (auto &entry : matchedMap)
        matchedPositions.push_back(entry.first);

    return matchedPositions.empty() ? STATUS_CODE_NOT_FOUND : STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KalmanSplittingAlgorithm::FollowRoute(const std::map<int, CaloHitList> &caloHitWireMap, KalmanFit &kalmanFit) const
{
    bool prevAmbiguous(true);

    // Continue fit
    for (auto &entry : caloHitWireMap)
    {
        if (entry.first <= kalmanFit.m_currentWireID)
            continue;
                  
        bool thisAmbiguous(entry.second.size() != 1);

        // If we've hit another ambiguity
        if (thisAmbiguous && !prevAmbiguous)
            return;

        // Make next kalman step
        kalmanFit.m_kalmanFilter.Predict();
       
        // Find matched cluster position
        CartesianPointVector matchedPositions;
        if (this->FindMatchedClusterPosition(entry.second, kalmanFit, matchedPositions) == STATUS_CODE_SUCCESS)
        {
            kalmanFit.SaveStep(entry.first);
            kalmanFit.AddPositionAndUpdate(matchedPositions.front());

            ////////////////////////////
            // CartesianVector bestPos(matchedPositions.front());
            // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &bestPos, "kalmanPos", VIOLET, 2));
            ////////////////////////////
        }           

        prevAmbiguous = thisAmbiguous;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanSplittingAlgorithm::FindBestSplitPosition(const Cluster *const pCluster, const KalmanFit &kalmanFit, 
    CartesianVector &splitPosition) const
{
    // Probe kalman fit
    const int nFitPoints(kalmanFit.m_positions.size());
    if (nFitPoints < ((2 * m_segmentWindows) + 1)) { return STATUS_CODE_NOT_FOUND; };

    bool splitFound(false);
    float maxDeviation(std::numeric_limits<float>::min());

    for (unsigned int i = (m_segmentWindows + 1); i < static_cast<unsigned int>(nFitPoints - m_segmentWindows); ++i)
    {
        // Get before/after segments
        CartesianPointVector beforeSeg, afterSeg;
        for (int j = 1; j <= m_segmentWindows; j++)
        {
            beforeSeg.push_back(kalmanFit.m_directions.at(i - j)); 
            afterSeg.push_back(kalmanFit.m_directions.at(i + j));
        }

        // Get median direction from before/after region
        const CartesianVector beforeMedian(this->GetMedian(beforeSeg));
        const CartesianVector afterMedian(this->GetMedian(afterSeg));

        // Get splitting metrucs
        const float eitherSide(beforeMedian.GetOpeningAngle(afterMedian) / 3.14 * 180.f);
        const float beforeSigma(this->GetSTD(beforeSeg)); // make sure that before segment looks track-like...

        // Measure angular deviation
        if ((eitherSide > m_minDeviation) && (beforeSigma < m_maxSpread))
        { 
            if (eitherSide > maxDeviation)
            {
                splitFound = true;
                maxDeviation = eitherSide;
                splitPosition = kalmanFit.m_positions.at(i);
            }
        }
   } 

    if (splitFound)
    {
        PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &splitPosition, "splitPosition", BLUE, 2));
        // PANDORA_MONITORING_API(ViewEvent(this->GetPandora()));
    }

    return splitFound ? STATUS_CODE_SUCCESS : STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector KalmanSplittingAlgorithm::GetMedian(CartesianPointVector &cartesianPointVector) const
{
    // Sort wrt x-axis theta
    std::sort(cartesianPointVector.begin(), cartesianPointVector.end(), [](const CartesianVector &lhs, const CartesianVector &rhs) 
        { 
            CartesianVector xAxis(1.f, 0.f, 0.f);
            float lhsAngle(xAxis.GetOpeningAngle(lhs));
            float rhsAngle(xAxis.GetOpeningAngle(rhs));

            if (lhs.GetY() < 0.f)
                lhsAngle += M_PI;

            if (rhs.GetY() < 0.f)
                rhsAngle += M_PI;

            return lhsAngle > rhsAngle;
        }
    );

    const int medianIndex = std::floor(cartesianPointVector.size() / 2);

    return cartesianPointVector.at(medianIndex);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float KalmanSplittingAlgorithm::GetSTD(const CartesianPointVector &directionVector) const
{
    CartesianVector xAxis(1.f, 0.f, 0.f);
    float mean(0.f), sigma(0.f);
    int count(0);

    for (const CartesianVector &direction : directionVector)
    {
        if (direction.GetMagnitude() < std::numeric_limits<float>::epsilon())
            continue;

        ++count;

        float openingAngle(xAxis.GetOpeningAngle(direction));

        if (direction.GetY() < 0.f)
            openingAngle += M_PI;

        mean += openingAngle;
    }

    mean /= count;

    for (const CartesianVector &direction : directionVector)
    {
        if (direction.GetMagnitude() < std::numeric_limits<float>::epsilon())
            continue;

        float openingAngle(xAxis.GetOpeningAngle(direction));

        if (direction.GetY() < 0.f)
            openingAngle += M_PI;

        sigma += std::pow(openingAngle - mean, 2);
    }

    sigma = std::sqrt(sigma / count);

    return sigma;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanSplittingAlgorithm::DivideCaloHits(const TwoDSlidingFitResult &slidingFitResult,
    const CartesianVector &splitPosition, CaloHitList &firstCaloHitList, CaloHitList &secondCaloHitList) const
{
    float rL(0.f), rT(0.f);
    slidingFitResult.GetLocalPosition(splitPosition, rL, rT);

    const Cluster *const pCluster(slidingFitResult.GetCluster());
    const OrderedCaloHitList &orderedCaloHitList(pCluster->GetOrderedCaloHitList());

    for (OrderedCaloHitList::const_iterator iter = orderedCaloHitList.begin(); iter != orderedCaloHitList.end(); ++iter)
    {
        for (CaloHitList::const_iterator hitIter = iter->second->begin(), hitIterEnd = iter->second->end(); hitIter != hitIterEnd; ++hitIter)
        {
            const CaloHit *const pCaloHit = *hitIter;

            float thisL(0.f), thisT(0.f);
            slidingFitResult.GetLocalPosition(pCaloHit->GetPositionVector(), thisL, thisT);

            if (thisL < rL)
            {
                firstCaloHitList.push_back(pCaloHit);
            }
            else
            {
                secondCaloHitList.push_back(pCaloHit);
            }
        }
    }

    if (firstCaloHitList.empty() || secondCaloHitList.empty())
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KalmanSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", m_minClusterLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "KalmanDelta", m_kalmanDelta));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "KalmanProcessVarCoeff", m_kalmanProcessVarCoeff));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "KalmanMeasurementVarCoeff", m_kalmanMeasurementVarCoeff));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinTransSeparation", m_minTransSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SegmentWindows", m_segmentWindows));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinDeviation", m_minDeviation));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSpread", m_maxSpread));

    return ClusterSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
