/**
 *  @file   larpandoracontent/LArHelpers/LArKalmanHelper.cc
 *
 *  @brief  Implementation of the kalman helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArHitWidthHelper.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/KalmanFit.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArKalmanHelper.h"
//#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "Managers/GeometryManager.h"
#include "Geometry/LArTPC.h"

#include <algorithm>
#include <cmath>
#include <limits>

using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------
    
KalmanFit LArKalmanHelper::PerformKalmanFit(const Pandora &pandora, const Cluster *const pCluster, const float minTransSep)
{
    // Get CaloHits
    CaloHitList caloHitList;
    LArClusterHelper::GetAllHits(pCluster, caloHitList);

    // Create wireID map to assess ambiguities
    const LArTPC *const pTPC(pandora.GetGeometry()->GetLArTPCMap().begin()->second);
    const HitType view(caloHitList.front()->GetHitType());
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
   
    std::map<int, CaloHitList> caloHitWireMap;
    for (const CaloHit *const pCaloHit : caloHitList)
        caloHitWireMap[std::floor(pCaloHit->GetPositionVector().GetZ() / pitch)].push_back(pCaloHit);

    // Initialise Kalman fit
    Eigen::VectorXd init(2);
    int startWireID(-1);

    if (LArKalmanHelper::GetFitSeed(caloHitWireMap, init, startWireID) != STATUS_CODE_SUCCESS)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);

    const float m_kalmanDelta(1.f), m_kalmanProcessVarCoeff(1.f), m_kalmanMeasurementVarCoeff(1.f);
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
        if (LArKalmanHelper::FindMatchedClusterPosition(entry.second, kalmanFit, minTransSep, matchedPositions) != STATUS_CODE_SUCCESS)
            continue;

        if (matchedPositions.size() == 1) // If unambiguous, just add in hit
        {
            kalmanFit.SaveStep(entry.first);
            kalmanFit.AddPositionAndUpdate(matchedPositions.front());

            ///////////////////////////////////
            // const CartesianVector bestPos(matchedPositions.front());
            // PANDORA_MONITORING_API(AddMarkerToVisualization(pandora, &bestPos, "kalmanPos", GREEN, 2));                        
            ///////////////////////////////////
        }
        else if (matchedPositions.size() > 1)  // If ambiguous, follow each route and pick the straightest
        {
            ///////////////////////////////////
            // PANDORA_MONITORING_API(ViewEvent(pandora));
            // for (const CartesianVector &seed : matchedPositions)
            // {
            //     std::cout << "Seed: " << seed << std::endl;
            //     PANDORA_MONITORING_API(AddMarkerToVisualization(pandora, &seed, "kalmanPos", RED, 2));
            // }
            // PANDORA_MONITORING_API(ViewEvent(pandora));
            ///////////////////////////////////

            KalmanFitVector pathways;

            for (const CartesianVector &seed : matchedPositions)
            {
                pathways.push_back(kalmanFit);
                pathways.back().SaveStep(entry.first);
                pathways.back().AddPositionAndUpdate(seed);

                // PANDORA_MONITORING_API(AddMarkerToVisualization(pandora, &seed, "kalmanPos", RED, 2));

                LArKalmanHelper::FollowRoute(caloHitWireMap, minTransSep, pathways.back());

                // PANDORA_MONITORING_API(ViewEvent(pandora));
            }

            // PANDORA_MONITORING_API(ViewEvent(pandora));

            // Find the best path
            int highestNHits(0);
            float straightest(std::numeric_limits<float>::max());

            for (unsigned int i = 0; i < pathways.size(); ++i)
            {
                const KalmanFit &pathwayFit(pathways.at(i));
                const int nHits(pathwayFit.m_positions.size());
                const float sigma(LArKalmanHelper::GetSTD(pathwayFit.m_directions));

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
            //     PANDORA_MONITORING_API(AddMarkerToVisualization(pandora, &frog, "CHOSEN", GREEN, 2));
            // }
            // PANDORA_MONITORING_API(ViewEvent(pandora));
            //////////////////////////////////
        }
    }

    //////////////////////////////////
    // PANDORA_MONITORING_API(ViewEvent(pandora));
    //////////////////////////////////
    
    return kalmanFit;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArKalmanHelper::GetFitSeed(const std::map<int, CaloHitList> &caloHitWireMap, Eigen::VectorXd &init, int &startWireID)
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

StatusCode LArKalmanHelper::FindMatchedClusterPosition(const CaloHitList &caloHitList, const KalmanFit &kalmanFit, 
    const float minTransSep, CartesianPointVector &matchedPositions)
{
    const CartesianVector kalmanPos(kalmanFit.m_kalmanFilter.GetPosition()(0), 0.f, kalmanFit.m_kalmanFilter.GetPosition()(1));
    const CartesianVector kalmanDir(kalmanFit.m_kalmanFilter.GetDirection()(0), 0.f, kalmanFit.m_kalmanFilter.GetDirection()(1));

    std::vector<std::pair<CartesianVector, float>> matchedMap;

    for (const CaloHit * const pCaloHit : caloHitList)
    {
        // if (pCaloHit->GetCellSize1() > 0.75f)
        //     continue;

        const CartesianVector hitPosition((pCaloHit->GetCellSize1() > 0.5f) ?
            LArHitWidthHelper::GetClosestPointToLine2D(kalmanPos, kalmanDir, pCaloHit) : pCaloHit->GetPositionVector()); 

        const float t(kalmanDir.GetCrossProduct(kalmanPos - hitPosition).GetMagnitude());

        if (t < minTransSep)
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

void LArKalmanHelper::FollowRoute(const std::map<int, CaloHitList> &caloHitWireMap, const float minTransSep, KalmanFit &kalmanFit)
{
    bool prevAmbiguous(true);

    // Continue fit
    for (auto &entry : caloHitWireMap)
    {
        if (entry.first <= kalmanFit.m_currentWireID)
            continue;
                 
        // Make next kalman step
        kalmanFit.m_kalmanFilter.Predict();
       
        // Find matched cluster position
        CartesianPointVector matchedPositions;
        if (LArKalmanHelper::FindMatchedClusterPosition(entry.second, kalmanFit, minTransSep, matchedPositions) != STATUS_CODE_SUCCESS)
            continue;

        bool thisAmbiguous(matchedPositions.size() != 1);

        // If we've hit another ambiguity
        if (thisAmbiguous && !prevAmbiguous)
            return;

        kalmanFit.SaveStep(entry.first);
        kalmanFit.AddPositionAndUpdate(matchedPositions.front());

        ////////////////////////////
        // CartesianVector bestPos(matchedPositions.front());
        // PANDORA_MONITORING_API(AddMarkerToVisualization(this->GetPandora(), &bestPos, "kalmanPos", VIOLET, 2));
        ////////////////////////////
                  
        prevAmbiguous = thisAmbiguous;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArKalmanHelper::GetSTD(const CartesianPointVector &directionVector)
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


} // namespace lar_content
