/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/CleanDLClusterSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the two dimensional sliding fit splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"
#include "larpandoracontent/LArUtility/KalmanFilter.h"
#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArTwoDReco/CleanDLClusterSplittingAlgorithm.h"

#include <torch/script.h>
#include <torch/torch.h>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

CleanDLClusterSplittingAlgorithm::CleanDLClusterSplittingAlgorithm() :
    m_nuVertexListName("NeutrinoVertices3D"),
    m_secVertexListName("SecondaryVertices3D"),
    m_pCaloHitList(nullptr),
    m_pClusterList(nullptr),
    m_pNuVertexList(nullptr),
    m_pSecVertexList(nullptr),
    m_minClusterHits(50),
    m_slidingWindow(20),
    m_lBinSize(0.5f),
    m_searchRegion1D(20.f),
    m_windowLength(48),
    m_isContaminatedThreshold(0.5f),
    m_isSplitThreshold(0.5f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CleanDLClusterSplittingAlgorithm::Run()
{
    if (this->GetLists() != STATUS_CODE_SUCCESS)
        return STATUS_CODE_SUCCESS;

    // For cluster features, we need a KDTree
    HitKDTree2D kdTree;
    HitKDNode2DList kdNode2DList;
    KDTreeBox kdTreeBox(fill_and_bound_2d_kd_tree(*m_pCaloHitList, kdNode2DList));
    kdTree.build(kdNode2DList, kdTreeBox);

    // Probe clusters
    ClusterVector internalClusterVector(m_pClusterList->begin(), m_pClusterList->end());
    std::sort(internalClusterVector.begin(), internalClusterVector.end(), LArClusterHelper::SortByNHits);

    for (const Cluster *const pCluster : internalClusterVector)
        this->ProcessCluster(pCluster, kdTree);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CleanDLClusterSplittingAlgorithm::GetLists()
{
    // Get 2D CaloHits - must find
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, m_pCaloHitList))

    if ((!m_pCaloHitList) || m_pCaloHitList->empty())
        return STATUS_CODE_NOT_FOUND;

    // Get 2D Clusters - must find
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_clusterListName, m_pClusterList))

    if ((!m_pClusterList) || m_pClusterList->empty())
        return STATUS_CODE_NOT_FOUND;

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, m_clusterListName))

    // Get NeutrinoVertex - must find
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_nuVertexListName, m_pNuVertexList));

    if (!m_pNuVertexList || (m_pNuVertexList->size() != 1))
        return STATUS_CODE_NOT_FOUND;

    // Get secondary vertices - okay if not found
    PandoraContentApi::GetList(*this, m_secVertexListName, m_pSecVertexList);

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CleanDLClusterSplittingAlgorithm::ProcessCluster(const Cluster *const pCluster, HitKDTree2D &kdTree)
{
    // Enough hits?
    CaloHitList clusterHits;
    LArClusterHelper::GetAllHits(pCluster, clusterHits);

    if (clusterHits.size() < m_minClusterHits)
        return;

    try
    {
        // Perform sliding fit
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        const TwoDSlidingFitResult clusterFit(pCluster, m_slidingWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));

        // Find pathway through the cluster
        ClusterPath clusterPath;
        this->FindPath(pCluster, clusterFit, clusterPath);
                      
        if (clusterPath.size() < m_windowLength)
            return;

        // Get features
        Features features;
        this->InitialiseFeatures(features);
        this->FillFeatures(clusterPath, clusterFit, features, kdTree);

        // Get split indices
        IntVector splitIndices;
        this->GetSplitIndices(features, splitIndices);

        // Remove any that are too close to the nu vertex (network picks up on hit sharing)
        this->FilterSplitIndices(clusterPath, hitType, splitIndices);

        if (splitIndices.empty())
            return;

        // Divide calo hits
        std::vector<CaloHitList> splitClusterHits(this->DivideCaloHits(clusterHits, clusterPath, clusterFit, splitIndices));

        if (splitClusterHits.size() < 2)
            return;

        // Split clusters
        this->SplitCluster(pCluster, splitClusterHits);
    }
    catch (...) {}
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CleanDLClusterSplittingAlgorithm::FindPath(const Cluster *const pCluster, const TwoDSlidingFitResult &clusterFit, ClusterPath &clusterPath) const
{
    // Get cluster hits
    CaloHitList clusterHits;
    LArClusterHelper::GetAllHits(pCluster, clusterHits);

    // Fit l decomposition
    for (const CaloHit *const pCaloHit : clusterHits)
    {
        float thisHitL(0.f), thisHitT(0.f);
        clusterFit.GetLocalPosition(pCaloHit->GetPositionVector(), thisHitL, thisHitT);
        const int lBinIndex(std::floor(thisHitL / m_lBinSize));

        ClusterHit clusterHit(pCaloHit, lBinIndex, thisHitT);

        const auto iter(std::find(clusterPath.begin(), clusterPath.end(), clusterHit));

        if (iter == clusterPath.end())
        {
            clusterPath.emplace_back(clusterHit);
        }
        else if (thisHitT < iter->m_t)
        {
            *iter = clusterHit;
        }
    }

    // Order low->high l
    std::sort(clusterPath.begin(), clusterPath.end(), [](const ClusterHit &lhs, const ClusterHit &rhs) { return lhs.m_l < rhs.m_l; });
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CleanDLClusterSplittingAlgorithm::InitialiseFeatures(Features &features) const
{
    features.insert(std::make_pair("Transverse", Feature((-0.01), 2.57, std::vector<float>())));
    features.insert(std::make_pair("Energy", Feature(0.49, 0.28, std::vector<float>())));
    features.insert(std::make_pair("HitWidth", Feature(0.56, 0.16, std::vector<float>())));
    features.insert(std::make_pair("Theta", Feature((-0.01), 0.11, std::vector<float>())));
    features.insert(std::make_pair("SecVertex", Feature(93.19, 103.39, std::vector<float>())));
    features.insert(std::make_pair("GapSep", Feature(171.78, 101.44, std::vector<float>())));
    features.insert(std::make_pair("EventHitSep", Feature(5.00, 5.83, std::vector<float>())));
    features.insert(std::make_pair("ClusterHitSep", Feature(0.54, 0.13, std::vector<float>())));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CleanDLClusterSplittingAlgorithm::FillFeatures(const ClusterPath &clusterPath, const TwoDSlidingFitResult &clusterFit, Features &features,
    HitKDTree2D &kdTree) const
{
    // Get path hits, and their total energy
    CaloHitList clusterPathHits; float totalEnergy(0.f);
    for (const auto &entry : clusterPath)
    {
        const CaloHit *const pPathCaloHit(entry.m_pHit);
        totalEnergy += pPathCaloHit->GetElectromagneticEnergy();
        clusterPathHits.push_back(pPathCaloHit);
    }

    // Initialise KalmanFilter
    KalmanFilter2D kalmanFilter2D(this->InitialiseKalmanFilter(clusterPath));

    // Fill features
    bool processedFirst(false); float cumulativeEnergy(0.f);
    for (unsigned int i = 0; i < clusterPath.size(); ++i)
    {
        const ClusterHit &clusterHit(clusterPath.at(i));

        if (!processedFirst)
        {
            features.at("Theta").m_sequence.push_back(-4.f);
            processedFirst = true;
        }
        else
        {
            // Make next step
            kalmanFilter2D.Predict();
            // Update feature
            Eigen::VectorXd eigenXd(2);
            const CartesianVector thisPosition(clusterHit.m_l, 0.f, clusterHit.m_t);
            eigenXd << thisPosition.GetX(), thisPosition.GetZ();
            kalmanFilter2D.Update(eigenXd);
            // Get scatter angle
            const CartesianVector newDirection(kalmanFilter2D.GetState()(2), 0.f, kalmanFilter2D.GetState()(3));
            float openingAngleL(-4.f);
            try
            {
                const float openingAngleT = CartesianVector(0.f, 0.f, 1.f).GetOpeningAngle(newDirection);
                openingAngleL = CartesianVector(1.f, 0.f, 0.f).GetOpeningAngle(newDirection);
                openingAngleL *= (openingAngleT > (M_PI * 0.5f)) ? (-1.f) : 1.f;
            }
            catch (...) {};
            features.at("Theta").m_sequence.push_back(openingAngleL);
        }

        const CaloHit *const pPathCaloHit(clusterHit.m_pHit);
        cumulativeEnergy += pPathCaloHit->GetElectromagneticEnergy();

        features.at("Transverse").m_sequence.push_back(clusterHit.m_t);
        features.at("Energy").m_sequence.push_back(cumulativeEnergy / totalEnergy);
        features.at("HitWidth").m_sequence.push_back(pPathCaloHit->GetCellSize1());
        features.at("SecVertex").m_sequence.push_back(
            this->GetDistanceToSecVertex(pPathCaloHit, pPathCaloHit->GetHitType()));
        features.at("EventHitSep").m_sequence.push_back(
            this->GetDistanceToEventHit(kdTree, pPathCaloHit, clusterPathHits));
        features.at("ClusterHitSep").m_sequence.push_back(i == (clusterPath.size() - 1) ? -1.f :
            this->GetDistanceToClusterHit(clusterHit, clusterPath.at(i + 1)));
        features.at("GapSep").m_sequence.push_back(i == 0 ? 1.f :
            this->GetDistanceToGap(clusterPath.at(i -1).m_pHit, pPathCaloHit));
    }

    // Smooth and normalise
    for (auto &entry : features)
    {
        if (entry.first == "GapSep")
            continue;

        entry.second.Smooth();
        entry.second.Normalise();
    }

    for (auto &entry : features)
    {
        std::cout << "----------" << std::endl;
        std::cout << entry.first << std::endl;
        for (int i =0; i < 5; ++i)
        {
            std::cout << entry.second.m_sequence.at(i) << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

KalmanFilter2D CleanDLClusterSplittingAlgorithm::InitialiseKalmanFilter(const ClusterPath &clusterPath) const
{
    // Kalman Config
    const LArTPC *const pTPC(this->GetPandora().GetGeometry()->GetLArTPCMap().begin()->second);
    const HitType view(clusterPath.begin()->m_pHit->GetHitType());
    const float pitch(view == TPC_VIEW_U ? pTPC->GetWirePitchU() : view == TPC_VIEW_V ? pTPC->GetWirePitchV() : pTPC->GetWirePitchW());
    const float kalmanDelta(1.f), kalmanProcessVarCoeff(1.f), kalmanMeasurementVarCoeff(1.f);
    const float processVariance{kalmanProcessVarCoeff * pitch * pitch};
    const float measurementVariance{kalmanMeasurementVarCoeff * pitch * pitch};

    // Initialise Kalman fit
    Eigen::VectorXd init(2);
    const float seedL(clusterPath.begin()->m_l), seedT(clusterPath.begin()->m_t);
    init << seedL, seedT;
    KalmanFilter2D kalmanFilter2D(kalmanDelta, processVariance, measurementVariance, init);

    return kalmanFilter2D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CleanDLClusterSplittingAlgorithm::GetAngle(const CartesianVector &position2D, const TwoDSlidingFitResult &clusterFit) const
{
    float thisL(0.f), thisT(0.f);
    clusterFit.GetLocalPosition(position2D, thisL, thisT);

    CartesianVector thisDirection(0.f, 0.f, 0.f);
    if (clusterFit.GetGlobalFitDirection(thisL, thisDirection) != STATUS_CODE_SUCCESS)
        return -4.f;

    const float openingAngleZ = CartesianVector(0.f, 0.f, 1.f).GetOpeningAngle(thisDirection);
    float openingAngleX = CartesianVector(1.f, 0.f, 0.f).GetOpeningAngle(thisDirection);
    openingAngleX *= (openingAngleZ > (M_PI * 0.5f)) ? (-1.f) : 1.f;

    return openingAngleX;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CleanDLClusterSplittingAlgorithm::GetDistanceToGap(const CaloHit *const pPrevHit, const CaloHit *const pCurrentHit) const
{
    if (!pPrevHit)
        return 1.f;

    const LArCaloHit *const pPrevLArHit(dynamic_cast<const LArCaloHit *>(pPrevHit));
    const LArCaloHit *const pCurrentLArHit(dynamic_cast<const LArCaloHit *>(pCurrentHit));

    if (!pPrevLArHit || !pCurrentLArHit)
        return 1.f;

    unsigned int prevTPCID(pPrevLArHit->GetLArTPCVolumeId());
    unsigned int currentTPCID(pCurrentLArHit->GetLArTPCVolumeId());

    if (prevTPCID != currentTPCID)
        return 0.f;

    unsigned int prevChildVolID(pPrevLArHit->GetDaughterVolumeId());
    unsigned int currentChildVolID(pCurrentLArHit->GetDaughterVolumeId());

    if (prevChildVolID != currentChildVolID)
        return 0.f;

    return 1.f;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CleanDLClusterSplittingAlgorithm::GetDistanceToEventHit(HitKDTree2D &kdTree, const CaloHit *const pCaloHit, const CaloHitList &clusterPathHits) const
{
     // Collect close hits
    HitKDNode2DList foundHits;
    KDTreeBox searchRegionHits(build_2d_kd_search_region(pCaloHit, m_searchRegion1D, m_searchRegion1D));
    kdTree.search(searchRegionHits, foundHits);
    // Filter
    bool found(false);
    float minDistSq(std::numeric_limits<float>::max());
    for (const auto &hit : foundHits)
    {
        const CaloHit *const pFoundHit(hit.data);

        if (std::find(clusterPathHits.begin(), clusterPathHits.end(), pFoundHit) != clusterPathHits.end())
            continue;

        found = true;
        minDistSq = std::min(minDistSq, (pCaloHit->GetPositionVector() - pFoundHit->GetPositionVector()).GetMagnitudeSquared());
    }

    return (found ? std::sqrt(minDistSq) : -1.f);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CleanDLClusterSplittingAlgorithm::GetDistanceToClusterHit(const ClusterHit &currentHit, 
    const ClusterHit &nextHit) const
{
    return (currentHit.m_pHit->GetPositionVector() - nextHit.m_pHit->GetPositionVector()).GetMagnitude();
}

//------------------------------------------------------------------------------------------------------------------------------------------

float CleanDLClusterSplittingAlgorithm::GetDistanceToSecVertex(const CaloHit *const pCaloHit, const HitType hitType) const
{
    if (!m_pSecVertexList)
        return -1.f;

    float bestSepSq(std::numeric_limits<float>::max());
    for (const Vertex *const pSecVertex : *m_pSecVertexList)
    {
        const CartesianVector secVtxPos(LArGeometryHelper::ProjectPosition(this->GetPandora(), pSecVertex->GetPosition(), hitType));
        bestSepSq = std::min(bestSepSq, (secVtxPos - pCaloHit->GetPositionVector()).GetMagnitudeSquared());
    }

    return std::sqrt(bestSepSq);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CleanDLClusterSplittingAlgorithm::GetSplitIndices(const Features &features, IntVector &splitIndices)
{
    // Find the window start indices
    IntVector windowStart;
    const int sequenceLength(features.at("Transverse").m_sequence.size());
    const int nWindows(std::floor(sequenceLength) / m_windowLength);

    for (int i = 0; i < nWindows; ++i)
        windowStart.push_back(m_windowLength * i);

    if (sequenceLength % m_windowLength != 0)
        windowStart.push_back(sequenceLength - m_windowLength);

    std::cout << "-------------" << std::endl;
    for (auto &entry : windowStart)
        std::cout << "windowStart: " << entry << std::endl;

    // Now get split indices and associated scores
    FloatVector splitScores;
    this->GetSplitIndices(features, windowStart, splitIndices, splitScores);

    // The model will often identify consecutive splitting positions around the truth, so identify one point
    this->FilterModelOutput(splitScores, splitIndices);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CleanDLClusterSplittingAlgorithm::GetSplitIndices(const Features &features, const IntVector &windowStart, IntVector &splitIndices, 
    FloatVector &splitScores)
{
    // Create network input
    LArDLHelper::TorchInput input; // (n_windows, sequence_length, n_features)
    LArDLHelper::InitialiseInput({static_cast<int>(windowStart.size()), m_windowLength, 8}, input); 

    for (unsigned int i = 0; i < windowStart.size(); ++i)
    {
        for (unsigned int j = 0; j < m_windowLength; ++j)
        {
            input[i][j][0] = features.at("Transverse").m_sequence.at(windowStart.at(i) + j);
            input[i][j][1] = features.at("Energy").m_sequence.at(windowStart.at(i) + j);
            input[i][j][2] = features.at("HitWidth").m_sequence.at(windowStart.at(i) + j);
            input[i][j][3] = features.at("Theta").m_sequence.at(windowStart.at(i) + j);
            input[i][j][4] = features.at("SecVertex").m_sequence.at(windowStart.at(i) + j);
            input[i][j][5] = features.at("GapSep").m_sequence.at(windowStart.at(i) + j);
            input[i][j][6] = features.at("EventHitSep").m_sequence.at(windowStart.at(i) + j);
            input[i][j][7] = features.at("ClusterHitSep").m_sequence.at(windowStart.at(i) + j);
        }
    }

    // Get model output
    LArDLHelper::TorchOutput windowOutput, splitPosOutput;
    LArDLHelper::Forward(m_windowModel, {input}, windowOutput);
    LArDLHelper::Forward(m_splitPosModel, {input}, splitPosOutput);
    torch::TensorAccessor<float, 2> windowOutputAccessor = windowOutput.accessor<float, 2>();
    torch::TensorAccessor<float, 3> splitPosOutputAccessor = splitPosOutput.accessor<float, 3>();

    for (unsigned int i = 0; i < windowStart.size(); ++i)
    {
        // Apply softmax
        float bkgProb(exp(windowOutputAccessor[i][0])), sigProb(exp(windowOutputAccessor[i][1])), shrProb(exp(windowOutputAccessor[i][2]));
        //float bkgProb_new = (bkgProb) / (bkgProb + sigProb + shrProb);
        float sigProb_new = (sigProb) / (bkgProb + sigProb + shrProb);
        //float shrProb_new = (shrProb) / (bkgProb + sigProb + shrProb);

        //std::cout << "--" << std::endl;
        //std::cout << "bkgProb: " << bkgProb_new << std::endl;
        //std::cout << "sigProb: " << sigProb_new << std::endl;
        //std::cout << "shrProb: " << shrProb_new << std::endl;
        //std::cout << "--" << std::endl;

        //std::cout << "isContamScore: " << sigProb << std::endl;

        // Is contaminated?
        if (sigProb_new > m_isContaminatedThreshold)
        {
            // Is there a split point?
            for (unsigned int j = 0; j < m_windowLength; ++j)
            {
                const int sequenceIndex(windowStart.at(i) + j);
                const int nextWindowStart((i+1) == windowStart.size() ? std::numeric_limits<int>::max() : windowStart.at(i+1));

                // make sure we only use the last window for the overlap region...
                if (sequenceIndex >= nextWindowStart)
                    continue;

                // Apply sigmoid
                float splitProb(1.f / (1.f + exp(-splitPosOutputAccessor[i][j][0])));

                std::cout << "splitProb: " << splitProb << std::endl;

                if (splitProb > m_isSplitThreshold)
                {
                    splitIndices.push_back(sequenceIndex);
                    splitScores.push_back(static_cast<float>(splitProb));
                }
            }
        }
    }

    // for (unsigned int i=0; i < splitIndices.size(); ++i)
    // {
    //     std::cout << "index: " << splitIndices.at(i) << std::endl;
    //     std::cout << "score: " << splitScores.at(i) << std::endl;
    // }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CleanDLClusterSplittingAlgorithm::FilterModelOutput(const FloatVector &splitScores, IntVector &splitIndices) const
{
    if (splitIndices.empty())
        return;

    IntVector splitIndices_temp(splitIndices);
    splitIndices.clear();

    int bestIndex = splitIndices_temp.front();
    float bestScore = splitScores.front();
    int previousIndex = splitIndices_temp.front();

    for (unsigned int i = 1; i < splitIndices_temp.size(); ++i)
    {
        int thisIndex = splitIndices_temp.at(i);
        float thisScore = splitScores.at(i);

        if (thisIndex == previousIndex + 1)
        {
            if (thisScore > bestScore)
            {
                bestScore = thisScore;
                bestIndex = thisIndex;
            }
        }
        else
        {
            splitIndices.push_back(bestIndex);

            bestIndex = thisIndex;
            bestScore = thisScore;
        }

        previousIndex = thisIndex;
    }

    // Push the last run
    splitIndices.push_back(bestIndex);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CleanDLClusterSplittingAlgorithm::FilterSplitIndices(const ClusterPath &clusterPath, const HitType hitType, IntVector &splitIndices) const
{
    IntVector filteredIndices;
    const CartesianVector &nuVertex3D(m_pNuVertexList->front()->GetPosition());
    const CartesianVector nuVertexPosition(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));

    for (const int splitIndex : splitIndices)
    {
        if (splitIndex >= static_cast<int>(clusterPath.size()))
        {
            std::cout << "AHHHHHHHH" << std::endl;
            throw;
        }

        const CartesianVector &position(clusterPath.at(splitIndex).m_pHit->GetPositionVector());
        
        if ((nuVertexPosition - position).GetMagnitude() > 5.f)
            filteredIndices.push_back(splitIndex);
    }

    splitIndices.swap(filteredIndices);
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::vector<CaloHitList> CleanDLClusterSplittingAlgorithm::DivideCaloHits(const CaloHitList &clusterHits, const ClusterPath &clusterPath, 
    const TwoDSlidingFitResult &clusterFit, const IntVector &splitIndices) const
{
    // Get l-coord of split index
    FloatVector lSplit;
    for (const int splitIndex : splitIndices)
    {
        CartesianVector position(clusterPath.at(splitIndex).m_pHit->GetPositionVector());
        float thisL(0.f), thisT(0.f);
        clusterFit.GetLocalPosition(position, thisL, thisT);
        lSplit.push_back(thisL);
     }

    // for (auto &entry : lSplit)
    // {
    //     std::cout << "lSplit: " << entry << std::endl;
    // }

    lSplit.insert(lSplit.begin(), std::numeric_limits<float>::lowest());
    lSplit.insert(lSplit.end(), std::numeric_limits<float>::max());
    const int nSplitPositions(lSplit.size() - 2);
    std::vector<CaloHitList> splitClusterHits((nSplitPositions + 1), CaloHitList());

    for (const CaloHit *const pCaloHit : clusterHits)
    {
        float thisL(0.f), thisT(0.f);
        clusterFit.GetLocalPosition(pCaloHit->GetPositionVector(), thisL, thisT);

        for (int i = 0; i <= nSplitPositions; ++i)
        {
            if ((thisL > lSplit.at(i)) && (thisL < lSplit.at(i + 1)))
            {
                splitClusterHits[i].push_back(pCaloHit);
                break;
            }
        }
    }

    // Remove empty lists
    splitClusterHits.erase(std::remove_if(splitClusterHits.begin(), splitClusterHits.end(),
        [](const CaloHitList& list){ return list.empty(); }), splitClusterHits.end());

    return splitClusterHits;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CleanDLClusterSplittingAlgorithm::SplitCluster(const Cluster *const pCluster, const std::vector<CaloHitList> &splitClusterHits) const
{
    // Begin cluster fragmentation operations
    const ClusterList clusterList(1, pCluster);
    std::string clusterListToSaveName, clusterListToDeleteName;

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::InitializeFragmentation(*this, clusterList, clusterListToDeleteName, clusterListToSaveName));

    for (const CaloHitList &caloHitList : splitClusterHits)
    {
        PandoraContentApi::Cluster::Parameters parameters;
        parameters.m_caloHitList = caloHitList;

        const Cluster *pNewCluster(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pNewCluster));
    }

    // End cluster fragmentation operations
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::EndFragmentation(*this, clusterListToSaveName, clusterListToDeleteName));
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CleanDLClusterSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, 
        XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, 
        XmlHelper::ReadValue(xmlHandle, "ClusterListName", m_clusterListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "SecVertexListName", m_secVertexListName));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "MinClusterHits", m_minClusterHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "SlidingWindow", m_slidingWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "LBinSize", m_lBinSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "SearchRegion1D", m_searchRegion1D));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "WindowLength", m_windowLength));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "IsContaminatedThreshold", m_isContaminatedThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "IsSplitThreshold", m_isSplitThreshold));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "WindowLength", m_windowLength));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "WindowModelName", m_windowModelName));
    m_windowModelName = LArFileHelper::FindFileInPath(m_windowModelName, "FW_SEARCH_PATH");
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_windowModelName, m_windowModel));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "SplitPosModelName", m_splitPosModelName));
    m_splitPosModelName = LArFileHelper::FindFileInPath(m_splitPosModelName, "FW_SEARCH_PATH");
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, LArDLHelper::LoadModel(m_splitPosModelName, m_splitPosModel));


    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
