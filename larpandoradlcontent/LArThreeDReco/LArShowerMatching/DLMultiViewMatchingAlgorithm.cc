/**
 *  @file   larpandoracontent/LArThreeDReco/LArThreeDBase/MultiViewMatchingAlgorithm.cc
 *
 *  @brief  Implementation of the n view matching algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArFileHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"


#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"
#include "larpandoradlcontent/LArHelpers/LArDLShowerHelper.h"

#include "larpandoradlcontent/LArThreeDReco/LArShowerMatching/DLMultiViewMatchingAlgorithm.h"


#include <chrono>


using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DLMultiViewMatchingAlgorithm::DLMultiViewMatchingAlgorithm() :
    m_trainingMode(false),
    m_trainingFileName("ShowerMatchingTraining.root"),
    m_trainingTreeName("trainingTree"),
    m_hitFeatureDim(14),
    m_polarRScaleFactor{1.f},
    m_cartesianXScaleFactor{1.f},
    m_cartesianZScaleFactor{1.f},
    m_matchThreshold(0.5f)    
{
}

//------------------------------------------------------------------------------------------------------------------------------------------    

DLMultiViewMatchingAlgorithm::WireID DLMultiViewMatchingAlgorithm::DecodeWireHash(const uint64_t id)
{
    int cryo(static_cast<unsigned int>((id >> 48) & 0xFFFF));
    int tpc(static_cast<unsigned int>((id >> 32) & 0xFFFF));
    int plane(static_cast<unsigned int>((id >> 16) & 0xFFFF));
    int wire(static_cast<unsigned int>(id & 0xFFFF));

    return WireID(cryo, tpc, plane, wire);
}    

//------------------------------------------------------------------------------------------------------------------------------------------    

StatusCode DLMultiViewMatchingAlgorithm::Run()
{
    // Get lists
    const ClusterList *pClusterListU(nullptr), *pClusterListV(nullptr), *pClusterListW(nullptr);
    if (this->GetList(pClusterListU, m_clusterListNameU) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;
    if (this->GetList(pClusterListV, m_clusterListNameV) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;
    if (this->GetList(pClusterListW, m_clusterListNameW) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;
    
    auto t0(std::chrono::high_resolution_clock::now());

    // Initialise
    ClusterExtentMap clusterExtentMap;
    this->PrepareClusters(pClusterListU, pClusterListV, pClusterListW, clusterExtentMap);
    
    auto t1(std::chrono::high_resolution_clock::now());
    
    // Get navigation maps
    this->FillNavigationMaps(clusterExtentMap);
    
    auto t2(std::chrono::high_resolution_clock::now());

    /////////////////////////////////////////////
    // View a navigation map
    // /////////////////////////////////////////////
    // for (auto &entry : navigationU)
    // {
    //     const ClusterList uVis({entry.first});
    //     const ClusterList projVis({entry.second});

    //     PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &uVis, "U cluster", RED);
    //     PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &projVis, "projected clusters", BLUE);
    //     PandoraMonitoringApi::ViewEvent(this->GetPandora());

    // }

    ClusterGroupVector clusterGroupVector;
    this->GetConnectedGroups(clusterGroupVector);

    /////////////////////////////////////////////        
    // View cluster groups
    /////////////////////////////////////////////    
    // for (const ClusterGroup &clusterGroup : clusterGroupVector)
    // {
    //     ClusterList visU(clusterGroup.m_clustersU);
    //     ClusterList visV(clusterGroup.m_clustersV);
    //     ClusterList visW(clusterGroup.m_clustersW);
    //     int nU(clusterGroup.m_clustersU.size()), nV(clusterGroup.m_clustersV.size()), nW(clusterGroup.m_clustersW.size());        

    //     std::cout << "BEFORE CLUSTER GROUP" << std::endl;
    //     std::cout << "nu: " << nU << ", nV: " << nV << ", nW: " << nW << std::endl;
    //     PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &visU, "U Clusters", RED);
    //     PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &visV, "V Clusters", BLUE);
    //     PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &visW, "W Clusters", BLACK);
    //     PandoraMonitoringApi::ViewEvent(this->GetPandora());         
    // }
    
    auto t3(std::chrono::high_resolution_clock::now());

    SimilarityMatrix globalSimMatrix;
    this->FillGlobalSimMatrix(clusterGroupVector, globalSimMatrix);

    auto t4(std::chrono::high_resolution_clock::now());     
    
    this->UpdateNavigationMaps(globalSimMatrix);

    /////////////////////////////////////////////
    // View sim matrix
    // ///////////////////////////////////////////// 
    // for (const auto &entry1 : globalSimMatrix)
    // {
    //     for (const auto &entry2 : entry1.second)
    //     {
    //         std::cout << "sim score: " << entry2.second << std::endl;
    //     }
    // }

    clusterGroupVector.clear();
    this->GetConnectedGroups(clusterGroupVector);    

    /////////////////////////////////////////////        
    // View cluster groups
    /////////////////////////////////////////////    
    for (const ClusterGroup &clusterGroup : clusterGroupVector)
    {
        ClusterList visU(clusterGroup.m_clustersU);
        ClusterList visV(clusterGroup.m_clustersV);
        ClusterList visW(clusterGroup.m_clustersW);
        int nU(clusterGroup.m_clustersU.size()), nV(clusterGroup.m_clustersV.size()), nW(clusterGroup.m_clustersW.size());        

        std::cout << "AFTER CLUSTER GROUP" << std::endl;
        std::cout << "nu: " << nU << ", nV: " << nV << ", nW: " << nW << std::endl;
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &visU, "U Clusters", RED);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &visV, "V Clusters", BLUE);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &visW, "W Clusters", BLACK);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());         
    }


    unsigned int repeatCounter(0);
    unsigned int m_nMaxRepeats(2);
    bool repeat(true);

    while (repeat && (repeatCounter < m_nMaxRepeats))
    {
        repeat = false;
        for (const auto &matchingTool : m_matchingToolVector)
        {
            const bool particlesMade(matchingTool->Run(this, globalSimMatrix));
            repeat = repeat ? repeat : particlesMade;
        }

        ++repeatCounter;
    }

    auto prepareT(std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0));
    auto navigationT(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1));
    auto connectionT(std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2));
    auto simMatrixconnectionT(std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3));      

    // std::cout << "prepare time: " << prepareT.count() << std::endl;
    // std::cout << "navigation map time: " << navigationT.count() << std::endl;
    // std::cout << "connected groups time: " << connectionT.count() << std::endl;
    // std::cout << "global sim matrix time: " << connectionT.count() << std::endl;   
    
    this->CleanUp();

    return STATUS_CODE_SUCCESS;
}


    

//------------------------------------------------------------------------------------------------------------------------------------------  

template <typename T>
StatusCode DLMultiViewMatchingAlgorithm::GetList(const T *&pList, const std::string listName)
{
    pList = nullptr;

    if (PandoraContentApi::GetList(*this, listName, pList) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_NOT_FOUND;

    if ((!pList) || pList->empty())
        return STATUS_CODE_NOT_FOUND;

    return STATUS_CODE_SUCCESS;
}


//------------------------------------------------------------------------------------------------------------------------------------------    

void DLMultiViewMatchingAlgorithm::PrepareClusters(const ClusterList *const pClusterListU, const ClusterList *const pClusterListV,
    const ClusterList *const pClusterListW, ClusterExtentMap &clusterExtentMap)
{
    for (const ClusterList *const pClusterList : {pClusterListU, pClusterListV, pClusterListW})
    {
        if (pClusterListU->empty())
            continue;

        const HitType hitType(LArClusterHelper::GetClusterHitType(pClusterList->front()));
        ClusterList &filteredClusters((hitType == TPC_VIEW_U) ? m_filteredU : (hitType == TPC_VIEW_V) ? m_filteredV : m_filteredW);
        
        for (const Cluster *const pCluster : *pClusterList)
        {
            if (!pCluster->IsAvailable())
                continue;
            
            if (pCluster->GetNCaloHits() < 5)
                continue;

            filteredClusters.emplace_back(pCluster);

            // Get drift-extent
            CartesianVector minCoord(0.f, 0.f, 0.f), maxCoord(0.f, 0.f, 0.f);
            LArClusterHelper::GetClusterBoundingBox(pCluster, minCoord, maxCoord);

            clusterExtentMap[pCluster] = std::pair(minCoord.GetX(), maxCoord.GetX());           
        }
    }    
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void DLMultiViewMatchingAlgorithm::FillNavigationMaps(const ClusterExtentMap &clusterExtentMap)
{
    std::map<HitType, std::vector<HitType>> projHitTypes({
            {TPC_VIEW_U, {TPC_VIEW_V, TPC_VIEW_W}},
            {TPC_VIEW_V, {TPC_VIEW_W}}
        });   

    for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V})
    {
        const ClusterList &clusters((hitType == TPC_VIEW_U) ? m_filteredU : (hitType == TPC_VIEW_V) ? m_filteredV : m_filteredW);
        NavigationMap &navigationMap((hitType == TPC_VIEW_U) ? m_navigationU : (hitType == TPC_VIEW_V) ? m_navigationV : m_navigationW);
        
        for (const HitType &projHitType : projHitTypes.at(hitType))
        {
            const ClusterList &projClusters((projHitType == TPC_VIEW_U) ? m_filteredU : (projHitType == TPC_VIEW_V) ? m_filteredV : m_filteredW);
            NavigationMap &projNavigationMap((projHitType == TPC_VIEW_U) ? m_navigationU : (projHitType == TPC_VIEW_V) ? m_navigationV : m_navigationW);

            for (const Cluster *const pCluster : clusters)
            {
                for (const Cluster *const pProjCluster : projClusters)
                {
                    if (this->DoClustersOverlap(pCluster, pProjCluster, clusterExtentMap))
                    {
                        navigationMap[pCluster].emplace_back(pProjCluster);
                        projNavigationMap[pProjCluster].emplace_back(pCluster);
                    }
                }
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------  
    
bool DLMultiViewMatchingAlgorithm::DoClustersOverlap(const Cluster *const pCluster1, const Cluster *const pCluster2, const ClusterExtentMap &clusterExtentMap)
{
    const std::pair<float, float> &driftExtent1(clusterExtentMap.at(pCluster1));
    const std::pair<float, float> &driftExtent2(clusterExtentMap.at(pCluster2));
    const float minX(std::max(driftExtent1.first, driftExtent2.first));
    const float maxX(std::min(driftExtent1.second, driftExtent2.second));            
    const float driftSpan(maxX - minX);
            
    if (driftSpan < std::numeric_limits<float>::epsilon())
        return false;

    return this->DoClustersOverlapInWire(pCluster1, pCluster2, minX, maxX);
}

//------------------------------------------------------------------------------------------------------------------------------------------      

bool DLMultiViewMatchingAlgorithm::DoClustersOverlapInWire(const Cluster *const pCluster1, const Cluster *const pCluster2, const float minX, const float maxX)
{
    CaloHitList caloHitList1, caloHitList2, filteredCaloHitList2;
    LArClusterHelper::GetAllHits(pCluster1, caloHitList1);
    LArClusterHelper::GetAllHits(pCluster2, caloHitList2);    
    for (const CaloHit *const pCaloHit2 : caloHitList2)
    {
        if ((pCaloHit2->GetPositionVector().GetX() < minX) || (pCaloHit2->GetPositionVector().GetX() > maxX))
            continue;

        filteredCaloHitList2.emplace_back(pCaloHit2);
    }

    unsigned int nSamplingPoints2(filteredCaloHitList2.size());

    // Search for overlap
    int overlapCount(0), nSamplingPoints1(0);
    for (const CaloHit *const pCaloHit1 : caloHitList1)
    {
        if ((pCaloHit1->GetPositionVector().GetX() < minX) || (pCaloHit1->GetPositionVector().GetX() > maxX))
            continue;

        ++nSamplingPoints1;
        
        for (const CaloHit *const pCaloHit2 : caloHitList2)
        {
            const LArCaloHit *const pLArHit1(dynamic_cast<const LArCaloHit *>(pCaloHit1));
            const WireID wireID1(this->DecodeWireHash(pLArHit1->GetWireHash()));            
            const WireID overlapMin11(this->DecodeWireHash(pLArHit1->GetOverlapMin1())); //plane 1
            const WireID overlapMax11(this->DecodeWireHash(pLArHit1->GetOverlapMax1()));
            const WireID overlapMin12(this->DecodeWireHash(pLArHit1->GetOverlapMin2())); //plane 2
            const WireID overlapMax12(this->DecodeWireHash(pLArHit1->GetOverlapMax2()));           
            const LArCaloHit *const pLArHit2(dynamic_cast<const LArCaloHit *>(pCaloHit2));            
            const WireID wireID2(this->DecodeWireHash(pLArHit2->GetWireHash()));

            // tpc should be the same
            if ((wireID1.m_cryostat != wireID2.m_cryostat) || (wireID1.m_tpc != wireID2.m_tpc))
                continue;

            if (wireID2.m_plane == overlapMin11.m_plane)
            {
                if ((overlapMin11.m_wire <= wireID2.m_wire) && (overlapMax11.m_wire >= wireID2.m_wire))
                {
                    overlapCount++;
                    break;
                }
            }
            else
            {
                if ((overlapMin12.m_wire <= wireID2.m_wire) && (overlapMax12.m_wire >= wireID2.m_wire))
                {
                    overlapCount++;
                    break;
                }
            }
        }
    }

    const float overlapFraction1(nSamplingPoints1 == 0 ? 0.f : float(overlapCount) / float(nSamplingPoints1));
    const float overlapFraction2(nSamplingPoints2 == 0 ? 0.f : float(overlapCount) / float(nSamplingPoints2));
    return ((overlapFraction1 > 0.8f) && (overlapFraction2 > 0.8f));
}


//------------------------------------------------------------------------------------------------------------------------------------------    

void DLMultiViewMatchingAlgorithm::GetConnectedGroups(ClusterGroupVector &clusterGroupVector)
{
    ClusterList usedU, usedV, usedW;
    
    for (const ClusterList &clusterList : {m_filteredU, m_filteredV, m_filteredW})
    {
        for (const Cluster *const pCluster : clusterList)
        {
            ClusterGroup clusterGroup;
            this->GetConnectedGroup(pCluster, clusterGroup, usedU, usedV, usedW);

            int nU(clusterGroup.m_clustersU.size()), nV(clusterGroup.m_clustersV.size()), nW(clusterGroup.m_clustersW.size());
            
            if ((nU + nV + nW) <= 1)
                continue;

            int n0(0);
            for (int nView : {nU, nV, nW})
                if (nView == 0)
                    ++n0;

            if (n0 >=2)
                std::cout << "THIS IS BAD!" << std::endl;

            clusterGroupVector.push_back(clusterGroup);
           
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void DLMultiViewMatchingAlgorithm::GetConnectedGroup(const Cluster *const pCluster, ClusterGroup &clusterGroup, ClusterList &usedU, ClusterList &usedV, ClusterList &usedW)
{
    if (!pCluster->IsAvailable())
        return;
    
    const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
    ClusterList &used((hitType == TPC_VIEW_U) ? usedU : (hitType == TPC_VIEW_V) ? usedV : usedW);
        
    if (std::find(used.begin(), used.end(), pCluster) != used.end())
        return;

    // Add cluster
    ClusterList &clusterGroupList((hitType == TPC_VIEW_U) ? clusterGroup.m_clustersU : (hitType == TPC_VIEW_V) ? clusterGroup.m_clustersV : clusterGroup.m_clustersW);
    clusterGroupList.emplace_back(pCluster);
    used.emplace_back(pCluster);
    
    // Now find its connections and add those
    const NavigationMap &navigation((hitType == TPC_VIEW_U) ? m_navigationU : (hitType == TPC_VIEW_V) ? m_navigationV : m_navigationW);
    const auto navIter(navigation.find(pCluster));

    if (navIter == navigation.end())
        return;
    
    for (const Cluster *const pNavCluster : navIter->second)
        this->GetConnectedGroup(pNavCluster, clusterGroup, usedU, usedV, usedW);
}

//------------------------------------------------------------------------------------------------------------------------------------------     
    
void DLMultiViewMatchingAlgorithm::FillGlobalSimMatrix(const ClusterGroupVector &clusterGroupVector, SimilarityMatrix &globalSimMatrix)
{
    // Get nu vertex projections
    const VertexList *pVertexList(nullptr);
    if (this->GetList(pVertexList, m_nuVertexListName) != STATUS_CODE_SUCCESS)
        return;
    const Vertex *const pVertex(pVertexList->front());
    PANDORA_THROW_IF(STATUS_CODE_INVALID_PARAMETER, pVertex->GetVertexType() != VERTEX_3D);
    const std::map<HitType, CartesianVector> viewToVtxPos(
           {{TPC_VIEW_U, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_U)},
            {TPC_VIEW_V, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_V)},
            {TPC_VIEW_W, LArGeometryHelper::ProjectPosition(this->GetPandora(), pVertex->GetPosition(), TPC_VIEW_W)}});

    // Get detector x-gaps
    std::set<float> detXGaps;
    LArGeometryHelper::GetDetectorXGaps(this->GetPandora(), detXGaps);

    // Walk through each connected group
    for (const ClusterGroup &clusterGroup : clusterGroupVector)
    {
        this->PredictClusterSimilarityMatrix(clusterGroup.m_clustersU, clusterGroup.m_clustersV, clusterGroup.m_clustersW,
                                             viewToVtxPos, detXGaps, globalSimMatrix);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLMultiViewMatchingAlgorithm::PredictClusterSimilarityMatrix(const ClusterList &clusterListU, const ClusterList &clusterListV, const ClusterList &clusterListW, 
    const std::map<HitType, CartesianVector> &viewToVtxPos, const std::set<float> &detXGaps, SimilarityMatrix &clusterSimMat)
{
    // First sort the lists
    ClusterVector clusterVecU(clusterListU.begin(), clusterListU.end());
    std::sort(clusterVecU.begin(), clusterVecU.end(), LArClusterHelper::SortByNHits);
    ClusterVector clusterVecV(clusterListV.begin(), clusterListV.end());
    std::sort(clusterVecV.begin(), clusterVecV.end(), LArClusterHelper::SortByNHits);
    ClusterVector clusterVecW(clusterListW.begin(), clusterListW.end());
    std::sort(clusterVecW.begin(), clusterVecW.end(), LArClusterHelper::SortByNHits);
    int nClusters(clusterVecU.size() + clusterVecV.size() + clusterVecW.size());
    ClusterVector clusterVec;
    clusterVec.insert(clusterVec.end(), clusterVecU.begin(), clusterVecU.end());
    clusterVec.insert(clusterVec.end(), clusterVecV.begin(), clusterVecV.end());
    clusterVec.insert(clusterVec.end(), clusterVecW.begin(), clusterVecW.end());

    // Get our similarity matrix
    torch::InferenceMode guard;
    std::vector<torch::Tensor> tensorEncodedClusters;

    for (const Cluster *const pCluster : clusterVec)
    {
        CaloHitList clusterHits;
        LArClusterHelper::GetAllHits(pCluster, clusterHits);

        const HitType view(LArClusterHelper::GetClusterHitType(pCluster));

        std::vector<LArDLShowerHelper::HitFeatures> clusterFeatures;
        for (const CaloHit *const pCaloHit : clusterHits)
        {
            LArDLShowerHelper::HitFeatures hitFeatures;
            LArDLShowerHelper::CalculateHitFeatures(pCaloHit, detXGaps, viewToVtxPos.at(view), hitFeatures);
            clusterFeatures.emplace_back(hitFeatures);
        }

        torch::Tensor tensorCluster;
        this->MakeClusterTensor(clusterFeatures, view, nClusters, tensorCluster);
        torch::Tensor tensorEncodedCluster{m_modelEncoder.forward({tensorCluster}).toTensor()};
        tensorEncodedClusters.emplace_back(tensorEncodedCluster);
    }

    torch::Tensor tensorEncodedEvent{torch::cat(tensorEncodedClusters, 1)};
    tensorEncodedClusters.clear(); // Free memory

    torch::Tensor tensorAttnEvent{m_modelAttn.forward({tensorEncodedEvent}).toTensor()};
    tensorEncodedEvent = torch::Tensor(); // Free memory

    torch::Tensor tensorSimMat{m_modelSim.forward({tensorAttnEvent}).toTensor()};
    tensorAttnEvent = torch::Tensor();

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, this->PopulateClusterSimilarityMatrix(tensorSimMat, clusterVec, clusterSimMat));

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLMultiViewMatchingAlgorithm::MakeClusterTensor(const std::vector<LArDLShowerHelper::HitFeatures> &clusterFeatures, const HitType view,
    const int nClusters, torch::Tensor &tensorCluster) const
{
    // Global features
    const int nHits{static_cast<int>(clusterFeatures.size())};
    const float nHitsFeat{std::log(static_cast<float>(nHits))};
    const float nClustersFeat{std::log(static_cast<float>(nClusters))};
    const float zWidthFeat{LArGeometryHelper::GetWirePitch(this->GetPandora(), view) * m_cartesianZScaleFactor};

    // Fill tensor
    tensorCluster = torch::zeros({1, nHits, m_hitFeatureDim});
    auto accessor = tensorCluster.accessor<float, 3>();
    for (int i = 0; i < nHits; i++)
    {
        const LArDLShowerHelper::HitFeatures hitFeatures{clusterFeatures.at(i)};
        accessor[0][i][0] = hitFeatures.m_rRel * m_polarRScaleFactor;
        accessor[0][i][1] = hitFeatures.m_cosThetaRel;
        accessor[0][i][2] = hitFeatures.m_sinThetaRel;
        accessor[0][i][3] = hitFeatures.m_xRel * m_cartesianXScaleFactor;
        accessor[0][i][4] = hitFeatures.m_zRel * m_cartesianZScaleFactor;
        accessor[0][i][5] = hitFeatures.m_xWidth * m_cartesianXScaleFactor;
        accessor[0][i][6] = zWidthFeat;
        accessor[0][i][7] = hitFeatures.m_distToXGap * m_cartesianXScaleFactor;
        accessor[0][i][8] = std::log(hitFeatures.m_energy);
        accessor[0][i][9] = view == TPC_VIEW_U ? 1.f : 0.f;
        accessor[0][i][10] = view == TPC_VIEW_V ? 1.f : 0.f;
        accessor[0][i][11] = view == TPC_VIEW_W ? 1.f : 0.f;
        accessor[0][i][12] = nHitsFeat;
        accessor[0][i][13] = nClustersFeat;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLMultiViewMatchingAlgorithm::PopulateClusterSimilarityMatrix(const torch::Tensor &tensorSimMat, const ClusterVector &clusterVector,
    SimilarityMatrix &clusterSimMat) const
{
    // Check predicted sim matrics
    //PANDORA_RETURN_IF(STATUS_CODE_NOT_ALLOWED, !clusterSimMat.empty());
    const int64_t nClusters{static_cast<int64_t>(clusterVector.size())};
    PANDORA_RETURN_IF(STATUS_CODE_NOT_ALLOWED, tensorSimMat.dim() != 3 || nClusters != tensorSimMat.size(-1) || nClusters != tensorSimMat.size(-2));
    
    // Populate SimilarityMatrix
    auto accessor = tensorSimMat.accessor<float, 3>();

    auto iterI{clusterVector.begin()};
    for (int i = 0; i < nClusters; i++, iterI++)
    {
        auto iterJ{clusterVector.begin()};
        for (int j = 0; j < nClusters; j++, iterJ++)
        {
            clusterSimMat[*iterI][*iterJ] = accessor[0][i][j];
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------    

void DLMultiViewMatchingAlgorithm::UpdateNavigationMaps(const SimilarityMatrix &globalSimMatrix)
{
    for (const HitType &hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        NavigationMap &navigationMap((hitType == TPC_VIEW_U) ? m_navigationU : (hitType == TPC_VIEW_V) ? m_navigationV : m_navigationW);

        for (auto &[pCluster, clusterList] : navigationMap)
        {
            for (auto iter = clusterList.begin(); iter != clusterList.end(); )
            {
                auto simIter1(globalSimMatrix.find(pCluster));
                if (simIter1 == globalSimMatrix.end()) {++iter; continue;}
                auto simIter2(simIter1->second.find(*iter));
                if (simIter2 == simIter1->second.end()) {++iter; continue;}                

                if (simIter2->second < m_matchThreshold)
                    iter = clusterList.erase(iter);
                else
                    ++iter;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLMultiViewMatchingAlgorithm::CreatePfo(const ClusterList &clusters)
{
    const PfoList *pPfoList(nullptr);
    std::string pfoListName;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::CreateTemporaryListAndSetCurrent(*this, pPfoList, pfoListName));

    PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
    pfoParameters.m_particleId = E_MINUS;
    pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
    pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
    pfoParameters.m_energy = 0.f;
    pfoParameters.m_momentum = CartesianVector(0.f, 0.f, 0.f);
    pfoParameters.m_clusterList = clusters;
    
    const ParticleFlowObject *pPfo(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));

    if (!pPfoList->empty())
    {
        std::string m_outputPfoListName("ShowerParticles3D");
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Pfo>(*this, m_outputPfoListName));
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Pfo>(*this, m_outputPfoListName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void DLMultiViewMatchingAlgorithm::DeleteCluster(const Cluster *const pClusterToRemove)
{
    // Remove from filtered lists
    const HitType hitType(LArClusterHelper::GetClusterHitType(pClusterToRemove));
    ClusterList &filteredClusters((hitType == TPC_VIEW_U) ? m_filteredU : (hitType == TPC_VIEW_V) ? m_filteredV : m_filteredW);

    auto iter(std::find(filteredClusters.begin(), filteredClusters.end(), pClusterToRemove));
    if (iter != filteredClusters.end()) { filteredClusters.erase(iter); }

    // Remove from navigation maps
    NavigationMap &navigation((hitType == TPC_VIEW_U) ? m_navigationU : (hitType == TPC_VIEW_V) ? m_navigationV : m_navigationW);
    auto iterNav(navigation.find(pClusterToRemove));

    if (iterNav != navigation.end())
    {
        for (const Cluster *const pAssocCluster : iterNav->second)
        {
            const HitType assocHitType(LArClusterHelper::GetClusterHitType(pAssocCluster));
            NavigationMap &assocNavigation((assocHitType == TPC_VIEW_U) ? m_navigationU : (assocHitType == TPC_VIEW_V) ? m_navigationV : m_navigationW);

            auto iterAssoc(assocNavigation.find(pAssocCluster));
            if (iterAssoc == assocNavigation.end()) { continue; }

            auto iterRemove(std::find(iterAssoc->second.begin(), iterAssoc->second.end(), pClusterToRemove));
            if (iterRemove != iterAssoc->second.end()) { iterAssoc->second.erase(iterRemove); }
        }

        navigation.erase(iterNav);
    }
}
    
//------------------------------------------------------------------------------------------------------------------------------------------    

// void DLMultiViewMatchingAlgorithm::CollectConnectedGroup(const Cluster *const pKeyClusterU, const NavigationMap &navigationMapUV,
//     const NavigationMap &navigationMapUW, const NavigationMap &navigationMapVW)
// {
//     ClusterList clusterGroup({pKeyClusterU});

//     // Go through V
//     if (navigationMapUV.find(pClusterU) != navigationMapUV.end())
//     {
//         for (const Cluster *const pClusterV : navigationMap.at()



// }


//------------------------------------------------------------------------------------------------------------------------------------------    

void DLMultiViewMatchingAlgorithm::CleanUp()
{
    m_filteredU.clear();
    m_filteredV.clear();
    m_filteredW.clear();
    m_navigationU.clear();
    m_navigationV.clear();
    m_navigationW.clear();
}

//------------------------------------------------------------------------------------------------------------------------------------------
    
StatusCode DLMultiViewMatchingAlgorithm::ReadSettings([[maybe_unused]] const TiXmlHandle xmlHandle)
{
    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, 
        XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "MatchingTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        DLShowerMatchingTool *const pMatchingTool(dynamic_cast<DLShowerMatchingTool *>(*iter));

        if (!pMatchingTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_matchingToolVector.push_back(pMatchingTool);
    }
    
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "NuVertexListName", m_nuVertexListName));
    if (m_nuVertexListName.empty())
        m_nuVertexListName = "NeutrinoVertices3D";
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameU", m_clusterListNameU));
    if (m_clusterListNameU.empty())
        m_clusterListNameU = "ClustersU";
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameV", m_clusterListNameV));
    if (m_clusterListNameV.empty())
        m_clusterListNameV = "ClustersV";
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameW", m_clusterListNameW));
    if (m_clusterListNameW.empty())
        m_clusterListNameW = "ClustersW";

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "TrainingMode", m_trainingMode));

    if (m_trainingMode)
    {
        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
            XmlHelper::ReadValue(xmlHandle, "TrainingFileName", m_trainingFileName));

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
            XmlHelper::ReadValue(xmlHandle, "TrainingTreeName", m_trainingTreeName));
    }
    else
    {
        std::string modelEncoderName, modelAttnName, modelSimName;
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelEncoderFileName", modelEncoderName));
        modelEncoderName = LArFileHelper::FindFileInPath(modelEncoderName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelEncoderName, m_modelEncoder);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelAttnFileName", modelAttnName));
        modelAttnName = LArFileHelper::FindFileInPath(modelAttnName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelAttnName, m_modelAttn);
        PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ModelSimFileName", modelSimName));
        modelSimName = LArFileHelper::FindFileInPath(modelSimName, "FW_SEARCH_PATH");
        LArDLHelper::LoadModel(modelSimName, m_modelSim);

        PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
            XmlHelper::ReadValue(xmlHandle, "HitFeatureDim", m_hitFeatureDim));
    }

    const LArGeometryHelper::DetectorBoundaries detBounds{LArGeometryHelper::GetDetectorBoundaries(this->GetPandora())};
    double xLow{static_cast<double>(detBounds.m_xBoundaries.first)}, xHigh{static_cast<double>(detBounds.m_xBoundaries.second)};
    double yLow{static_cast<double>(detBounds.m_yBoundaries.first)}, yHigh{static_cast<double>(detBounds.m_yBoundaries.second)};
    double zLow{static_cast<double>(detBounds.m_zBoundaries.first)}, zHigh{static_cast<double>(detBounds.m_zBoundaries.second)};
    m_polarRScaleFactor = static_cast<float>(1. / std::sqrt(std::pow(xHigh - xLow, 2.) + std::pow(yHigh - yLow, 2.) + std::pow(zHigh - zLow, 2.)));
    m_cartesianXScaleFactor = static_cast<float>(1. / (xHigh - xLow));
    m_cartesianZScaleFactor = static_cast<float>(1. / (zHigh - zLow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MatchThreshold", m_matchThreshold));
    
    return STATUS_CODE_SUCCESS;

}   

//------------------------------------------------------------------------------------------------------------------------------------------

template StatusCode DLMultiViewMatchingAlgorithm::GetList(const ClusterList *&, const std::string);
template StatusCode DLMultiViewMatchingAlgorithm::GetList(const MCParticleList *&, const std::string);

} // namespace lar_dl_content
