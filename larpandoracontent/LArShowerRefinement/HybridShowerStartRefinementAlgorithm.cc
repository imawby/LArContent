/**
 *  @file   larpandoracontent/LArShowerRefinement/HybridShowerStartRefinementAlgorithm.cc
 *
 *  @brief  Implementation of the shower start refinement base algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMonitoringHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArThreeDSlidingFitResult.h"
#include "larpandoracontent/LArObjects/LArPointingCluster.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"
#include "larpandoracontent/LArShowerRefinement/HybridShowerStartRefinementAlgorithm.h"
#include "larpandoracontent/LArShowerRefinement/HybridShowerStartRefinementBaseTool.h"

using namespace pandora;

namespace lar_content
{

HybridShowerStartRefinementAlgorithm::HybridShowerStartRefinementAlgorithm() : 
    m_binSize(0.005),     
    m_minElectronCompleteness(0.33f),
    m_minElectronPurity(0.5f),
    m_minGammaCompleteness(0.33f),
    m_thresholdSignalGammaDisplacement(-3.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HybridShowerStartRefinementAlgorithm::Run()
{
    PfoVector pfoVector;
    this->FillPfoVector(pfoVector);

    CartesianVector nuVertexPosition(0.f, 0.f, 0.f);
    if (this->GetNeutrinoVertex(nuVertexPosition) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_SUCCESS;

    this->FillGammaHitMap();
    this->FillElectronHitMap();

    const CaloHitList *pCaloHitListU(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitListU", pCaloHitListU));

    const CaloHitList *pCaloHitListV(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitListV", pCaloHitListV));

    const CaloHitList *pCaloHitListW(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitListW", pCaloHitListW));

    // run tools
    for (const ParticleFlowObject *const pPfo : pfoVector)
    {
        for (HybridShowerStartRefinementBaseTool *const pShowerStartRefinementTool : m_algorithmToolVector)
        {
            if (std::find(m_deletedPfos.begin(), m_deletedPfos.end(), pPfo) != m_deletedPfos.end())
                continue;

            pShowerStartRefinementTool->Run(this, pPfo, nuVertexPosition, pCaloHitListU, pCaloHitListV, pCaloHitListW);
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HybridShowerStartRefinementAlgorithm::FillPfoVector(PfoVector &pfoVector)
{
    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        if (PandoraContentApi::GetList(*this, pfoListName, pPfoList) != STATUS_CODE_SUCCESS)
            continue;

        if (!pPfoList || pPfoList->empty())
        {
            std::cout << "HybridShowerStartRefinementAlgorithm: unable to find pfo list " << pfoListName << std::endl;
            continue;
        }

        pfoVector.insert(pfoVector.begin(), pPfoList->begin(), pPfoList->end());
    }

    // This ordering is important.
    std::sort(pfoVector.begin(), pfoVector.end(), LArPfoHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HybridShowerStartRefinementAlgorithm::GetNeutrinoVertex(CartesianVector &neutrinoVertex)
{
    const VertexList *pNuVertexList(nullptr);
    const StatusCode statusCode(PandoraContentApi::GetList(*this, m_neutrinoVertexListName, pNuVertexList));

    if (statusCode != STATUS_CODE_SUCCESS)
        return statusCode;

    if (!pNuVertexList || (pNuVertexList->size() != 1))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "HybridShowerStartRefinementAlgorithm: unable to find vertex list " << m_neutrinoVertexListName << " if it does exist, it may have more than one nu vertex" << std::endl;

        return STATUS_CODE_NOT_INITIALIZED;
    }

    neutrinoVertex = pNuVertexList->front()->GetPosition();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList HybridShowerStartRefinementAlgorithm::GetAllHitsOfType(const HitType hitType)
{
    CaloHitList viewHitList;

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "HybridShowerStartRefinementBaseTool: unable to find calo hit list " << "CaloHitList2D" << std::endl;

        return viewHitList;
    }

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (pCaloHit->GetHitType() == hitType)
            viewHitList.push_back(pCaloHit);
    }

    return viewHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CaloHitList HybridShowerStartRefinementAlgorithm::GetXIntervalHitsOfType(const ParticleFlowObject *const pShowerPfo, const HitType hitType)
{
    CaloHitList intervalHitList;

    ClusterList clustersU, clustersV, clustersW;
    LArPfoHelper::GetClusters(pShowerPfo, TPC_VIEW_U, clustersU); 
    LArPfoHelper::GetClusters(pShowerPfo, TPC_VIEW_V, clustersV); 
    LArPfoHelper::GetClusters(pShowerPfo, TPC_VIEW_W, clustersW); 

    CartesianVector uMin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    CartesianVector uMax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    CartesianVector vMin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    CartesianVector vMax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    CartesianVector wMin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
    CartesianVector wMax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

    if (!clustersU.empty())
        LArClusterHelper::GetClusterBoundingBox(clustersU.front(), uMin, uMax);

    if (!clustersV.empty())
        LArClusterHelper::GetClusterBoundingBox(clustersV.front(), vMin, vMax);

    if (!clustersW.empty())
        LArClusterHelper::GetClusterBoundingBox(clustersW.front(), wMin, wMax);

    float xMin(std::min(std::min(uMin.GetX(), vMin.GetX()), wMin.GetX()));
    float xMax(std::max(std::max(uMax.GetX(), vMax.GetX()), wMax.GetX()));

    CaloHitList viewHitList;

    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "HybridShowerStartRefinementBaseTool: unable to find calo hit list " << "CaloHitList2D" << std::endl;

        return intervalHitList;
    }

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (pCaloHit->GetHitType() != hitType)
            continue;

        const CartesianVector &hitPosition(pCaloHit->GetPositionVector());

        if ((hitPosition.GetX() > xMin) && (hitPosition.GetX() < xMax))
            intervalHitList.push_back(pCaloHit);
    }

    return intervalHitList;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HybridShowerStartRefinementAlgorithm::FillOwnershipMaps()
{
    m_hitToClusterMapU.clear(); m_hitToClusterMapV.clear(); m_hitToClusterMapW.clear();
    m_clusterToPfoMapU.clear(); m_clusterToPfoMapV.clear(); m_clusterToPfoMapW.clear();

    // First fill pfo maps
    PfoVector pfoVector;

    for (const std::string &pfoListName : m_pfoListNames)
    {
        const PfoList *pPfoList(nullptr);
        //PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, );
        PandoraContentApi::GetList(*this, pfoListName, pPfoList);

        if (!pPfoList || pPfoList->empty())
        {
            std::cout << "HybridShowerStartRefinementAlgorithm: unable to find pfo list " << pfoListName << std::endl;
            continue;
        }

        pfoVector.insert(pfoVector.begin(), pPfoList->begin(), pPfoList->end());
    }


    for (const ParticleFlowObject *const pPfo : pfoVector)
    {
        ClusterList twoDClusterList;
        LArPfoHelper::GetClusters(pPfo, TPC_VIEW_U, twoDClusterList);
        LArPfoHelper::GetClusters(pPfo, TPC_VIEW_V, twoDClusterList);
        LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, twoDClusterList);

        for (const Cluster *const pCluster : twoDClusterList)
        {
            CaloHitList caloHitList;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

            CaloHitList isolated(pCluster->GetIsolatedCaloHitList());

            for (const CaloHit *const pIsolated : isolated)
            {
                if (std::find(caloHitList.begin(), caloHitList.end(), pIsolated) == caloHitList.end())
                    caloHitList.push_back(pIsolated);
            }

            const HitType hitType(caloHitList.front()->GetHitType());
            HitToClusterMap &hitToClusterMap(hitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);
            ClusterToPfoMap &clusterToPfoMap(hitType == TPC_VIEW_U ? m_clusterToPfoMapU : hitType == TPC_VIEW_V ? m_clusterToPfoMapV : m_clusterToPfoMapW);

            for (const CaloHit *const pCaloHit : caloHitList)
                hitToClusterMap[pCaloHit] = pCluster;

            clusterToPfoMap[pCluster] = pPfo;
        }
    }

    // Now fill cluster maps
    StringVector clusterListNames;
    clusterListNames.push_back("ClustersU");
    clusterListNames.push_back("ClustersV");
    clusterListNames.push_back("ClustersW");

    ClusterList clusterList;
    for (const std::string &clusterListName : clusterListNames)
    {
        const ClusterList *pClusterList(nullptr);
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, clusterListName, pClusterList));

        if (!pClusterList || pClusterList->empty())
        {
            std::cout << "HybridShowerStartRefinementAlgorithm: No cluster list found, returning..." << std::endl;
            throw;
        }

        clusterList.insert(clusterList.end(), pClusterList->begin(), pClusterList->end());
    }

    for (const Cluster *const pCluster : clusterList)
    {
        const HitType hitType(LArClusterHelper::GetClusterHitType(pCluster));
        HitToClusterMap &hitToClusterMap(hitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);
        ClusterToPfoMap &clusterToPfoMap(hitType == TPC_VIEW_U ? m_clusterToPfoMapU : hitType == TPC_VIEW_V ? m_clusterToPfoMapV : m_clusterToPfoMapW);

        if (clusterToPfoMap.find(pCluster) != clusterToPfoMap.end())
            continue;

        if (!pCluster->IsAvailable())
        {
            std::cout << "CLUSTER IS NOT AVAILABLE ISOBE" << std::endl;
            throw;
        }

        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        CaloHitList isolated(pCluster->GetIsolatedCaloHitList());

        for (const CaloHit *const pIsolated : isolated)
        {
            if (std::find(caloHitList.begin(), caloHitList.end(), pIsolated) == caloHitList.end())
                caloHitList.push_back(pIsolated);
        }

        for (const CaloHit *const pCaloHit : caloHitList)
            hitToClusterMap[pCaloHit] = pCluster;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HybridShowerStartRefinementAlgorithm::FillGammaHitMap()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "HybridShowerStartRefinementBaseTool: unable to find calo hit list " << "CaloHitList2D" << std::endl;

        return;
    }

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        try
        {
            const MCParticle *pMCParticle(MCParticleHelper::GetMainMCParticle(pCaloHit));
            
            if (pMCParticle->GetParticleId() == 22)
                m_gammaHitMap[pMCParticle].push_back(pCaloHit);
        }
        catch (...)
        {
            continue;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HybridShowerStartRefinementAlgorithm::FillElectronHitMap()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "CaloHitList2D", pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "HybridShowerStartRefinementBaseTool: unable to find calo hit list " << "CaloHitList2D" << std::endl;

        return;
    }

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        MCParticleVector contributingMCParticleVector;
        const MCParticleWeightMap &weightMap(pCaloHit->GetMCParticleWeightMap());

        for (const auto &mapEntry : weightMap)
            contributingMCParticleVector.push_back(mapEntry.first);

        std::sort(contributingMCParticleVector.begin(), contributingMCParticleVector.end(), PointerLessThan<MCParticle>());

        float highestWeight(0.f);
        const MCParticle *highestElectronContributor(nullptr);

        for (const MCParticle *const pMCParticle : contributingMCParticleVector)
        {
            const bool isLeadingElectron((std::abs(pMCParticle->GetParticleId()) == 11) && (LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle) == pMCParticle));

            if (isLeadingElectron)
            {
                const float weight(weightMap.at(pMCParticle));

                if (weight > highestWeight)
                {
                    highestWeight = weight;
                    highestElectronContributor = pMCParticle;
                }
            }
        }

        if (highestElectronContributor)
            m_electronHitMap[highestElectronContributor].push_back(pCaloHit);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HybridShowerStartRefinementAlgorithm::IsElectron(const ParticleFlowObject *const pPfo) const
{
    MCParticleVector mcElectronVector;

    for (auto &entry : m_electronHitMap)
        mcElectronVector.push_back(entry.first);

    CaloHitList pfoHitList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, pfoHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, pfoHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, pfoHitList);

    for (const MCParticle *const pMCElectron : mcElectronVector)
    {
        const CaloHitList &mcElectronHitList(m_electronHitMap.at(pMCElectron));
        const CaloHitList sharedHitList(LArMCParticleHelper::GetSharedHits(pfoHitList, mcElectronHitList));

        const float completeness(static_cast<float>(sharedHitList.size()) / static_cast<float>(mcElectronHitList.size()));
        const float purity(static_cast<float>(sharedHitList.size()) / static_cast<float>(pfoHitList.size()));

        if (completeness < m_minElectronCompleteness)
            continue;

        if (purity < m_minElectronPurity)
            continue;

        return true;
    }

    return false;
}

//------------------------------------------------------------------------------------------------------------------------------------------

// our signal here is actually gammas that have made a mistake by getting back to the nu vertex. i.e. their reco vertex is closer to the nu vertex than the truth says

bool HybridShowerStartRefinementAlgorithm::IsGamma(const ParticleFlowObject *const pPfo, const CartesianVector &nuVertexPosition) const
{
    const MCParticle *pMainMCParticle(LArMCParticleHelper::GetMainMCParticle(pPfo));
    int pdg(std::abs(pMainMCParticle->GetParticleId()));

    if (pdg == 11)
        return false;

    // but does it contain a large chunk of a gamma even if it is not the main owner?
    // find gamma with highest number of shared hits that has a completeness of >50 % in one view
    if (pdg != 22)
    {
        unsigned int highestSharedHits(0);

        MCParticleVector mcGammaVector;
        for (auto &entry : m_gammaHitMap)
        {
            if (entry.first->GetParticleId() == 22)
            {
                const CaloHitList &mcGammaHitList(m_gammaHitMap.at(entry.first));

                if (mcGammaHitList.size() > 100)
                    mcGammaVector.push_back(entry.first);
            }
        }

        std::sort(mcGammaVector.begin(), mcGammaVector.end(), LArMCParticleHelper::SortByMomentum);

        CaloHitList pfoHitList;
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, pfoHitList);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, pfoHitList);
        LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, pfoHitList);

        for (const MCParticle *const pMCGamma : mcGammaVector)
        {
            const CaloHitList &mcGammaHitList(m_gammaHitMap.at(pMCGamma));
            const CaloHitList sharedHitList(LArMCParticleHelper::GetSharedHits(pfoHitList, mcGammaHitList));

           // get each view completeness
            for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
            {
                const float completeness(static_cast<float>(LArMonitoringHelper::CountHitsByType(hitType, sharedHitList)) /
                    static_cast<float>(LArMonitoringHelper::CountHitsByType(hitType, mcGammaHitList)));

                if ((completeness > m_minGammaCompleteness) && (sharedHitList.size() > highestSharedHits))
                {
                    pMainMCParticle = pMCGamma;
                    highestSharedHits = sharedHitList.size();
                }
            }
        }
    }

    pdg = std::abs(pMainMCParticle->GetParticleId());

    if (pdg != 22)
        return false;

    CaloHitList caloHitList3D;
    LArPfoHelper::GetCaloHits(pPfo, TPC_3D, caloHitList3D);

    if (caloHitList3D.empty())
        return false;

    float closestDistance(std::numeric_limits<float>::max());
    CartesianVector showerVertex(0.f, 0.f, 0.f);

    for (const CaloHit *const pCaloHit : caloHitList3D)
    {
        float jam = (pCaloHit->GetPositionVector() - nuVertexPosition).GetMagnitude();
        if (jam < closestDistance)
        {
            closestDistance = jam;
            showerVertex = pCaloHit->GetPositionVector();
        }
    }

    const CartesianVector &mcVertex(pMainMCParticle->GetVertex());
    const float mcSeparation((mcVertex - nuVertexPosition).GetMagnitude());
    const float recoSeparation((showerVertex - nuVertexPosition).GetMagnitude());
    const float difference(recoSeparation - mcSeparation);

    //////////////////////////////////////
    /*
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &mcVertex, "mcVertex", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &showerVertex, "showerVertex", GREEN, 2);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    //////////////////////////////////////

    if (difference > m_thresholdSignalGammaDisplacement)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------


void HybridShowerStartRefinementAlgorithm::AddElectronPathway(const ParticleFlowObject *const pShowerPfo, const CaloHitList &pathwayHitList)
{
    // This is so incredibly lazy isobel
    this->FillOwnershipMaps();

    const HitType hitType(pathwayHitList.front()->GetHitType());

    ClusterList showerClusters2D, showerClusters3D;
    LArPfoHelper::GetClusters(pShowerPfo, hitType, showerClusters2D);
    LArPfoHelper::GetClusters(pShowerPfo, TPC_3D, showerClusters3D);

    const ThreeDSlidingFitResult showerSlidingFit(showerClusters3D.front(), 1000, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
    const CartesianVector &showerDirection3D(showerSlidingFit.GetGlobalMaxLayerDirection());

    std::map<const ParticleFlowObject*, int> showerHitCountMap;

    for (const CaloHit *const pPathwayHit : pathwayHitList)
    {
        const HitToClusterMap &hitToClusterMap(hitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);
        const ClusterToPfoMap &clusterToPfoMap(hitType == TPC_VIEW_U ? m_clusterToPfoMapU : hitType == TPC_VIEW_V ? m_clusterToPfoMapV : m_clusterToPfoMapW);

        if (hitToClusterMap.find(pPathwayHit) == hitToClusterMap.end())
            continue;

        if (clusterToPfoMap.find(hitToClusterMap.at(pPathwayHit)) == clusterToPfoMap.end())
            continue;

        const ParticleFlowObject *const pPathwayShower(clusterToPfoMap.at(hitToClusterMap.at(pPathwayHit)));

        if (LArPfoHelper::IsTrack(pPathwayShower))
            continue;

        if (pPathwayShower == pShowerPfo)
            continue;

        if (showerHitCountMap.find(pPathwayShower) == showerHitCountMap.end())
            showerHitCountMap[pPathwayShower] = 1;
        else
            showerHitCountMap[pPathwayShower] = showerHitCountMap[pPathwayShower] + 1;
    }

    PfoList significantShowersToMerge;

    for (const auto &entry : showerHitCountMap)
    {
        float contaminationRatio(static_cast<float>(entry.second) / static_cast<float>(pathwayHitList.size()));

        if (contaminationRatio < 0.3f)
            continue;

        ClusterList pathwayClusters3D;
        LArPfoHelper::GetClusters(entry.first, TPC_3D, pathwayClusters3D);

        if (pathwayClusters3D.size() == 0)
            continue;

        try
        {
            const ThreeDSlidingFitResult pathwaySlidingFit(pathwayClusters3D.front(), 1000, LArGeometryHelper::GetWireZPitch(this->GetPandora()));
            const CartesianVector &pathwayDirection(pathwaySlidingFit.GetGlobalMaxLayerDirection());

            const float openingAngle(pathwayDirection.GetOpeningAngle(showerDirection3D) * 180 / M_PI);

            if (openingAngle < 5.f)
                significantShowersToMerge.push_back(entry.first);
        }
        catch (...)
        {
        }
    }

    // Add in hits first, then deal with merges
    for (const CaloHit *const pPathwayHit : pathwayHitList)
    {
        const HitToClusterMap &hitToClusterMap(hitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);
        const ClusterToPfoMap &clusterToPfoMap(hitType == TPC_VIEW_U ? m_clusterToPfoMapU : hitType == TPC_VIEW_V ? m_clusterToPfoMapV : m_clusterToPfoMapW);
        std::string clusterListName(hitType == TPC_VIEW_U ? "ClustersU" : hitType == TPC_VIEW_V ? "ClustersV" : "ClustersW");

        const Cluster *pParentCluster(nullptr);
        const ParticleFlowObject *pParentPfo(nullptr);

        if (hitToClusterMap.find(pPathwayHit) != hitToClusterMap.end())
        {
            pParentCluster = hitToClusterMap.at(pPathwayHit);

            if (clusterToPfoMap.find(pParentCluster) != clusterToPfoMap.end())
            {
                pParentPfo = clusterToPfoMap.at(pParentCluster);

                if (pParentPfo == pShowerPfo)
                    continue;

                if (std::find(significantShowersToMerge.begin(), significantShowersToMerge.end(), pParentPfo) != significantShowersToMerge.end())
                    continue;
            }
        }

        if (pParentCluster)
        {
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));

            CaloHitList clusterNormalHitList; const CaloHitList clusterIsolatedHitList(pParentCluster->GetIsolatedCaloHitList());
            pParentCluster->GetOrderedCaloHitList().FillCaloHitList(clusterNormalHitList);

            const bool isIsolated(std::find(clusterIsolatedHitList.begin(), clusterIsolatedHitList.end(), pPathwayHit) != clusterIsolatedHitList.end());

            if (!isIsolated && (clusterNormalHitList.size() == 1) && !(clusterIsolatedHitList.empty()))
            {
                const HitType isolatedHitType(LArClusterHelper::GetClusterHitType(pParentCluster));
                HitToClusterMap &isolatedHitToClusterMap(isolatedHitType == TPC_VIEW_U ? m_hitToClusterMapU : hitType == TPC_VIEW_V ? m_hitToClusterMapV : m_hitToClusterMapW);

                for (const CaloHit * const pIsolatedHit : clusterIsolatedHitList)
                {
                    isolatedHitToClusterMap.erase(pIsolatedHit);
                    const StatusCode isolatedStatusCode(PandoraContentApi::RemoveIsolatedFromCluster(*this, pParentCluster, pIsolatedHit));

                    if (isolatedStatusCode != STATUS_CODE_SUCCESS)
                    {
                        std::cout << "ISOBEL CANNOT REMOVE ISOLATED HIT?" << std::endl;
                        throw;
                    }
                }
            }

            const StatusCode statusCodeCluster(isIsolated ? PandoraContentApi::RemoveIsolatedFromCluster(*this, pParentCluster, pPathwayHit) : 
                PandoraContentApi::RemoveFromCluster(*this, pParentCluster, pPathwayHit));

            if (statusCodeCluster != STATUS_CODE_SUCCESS)
            {
                if (statusCodeCluster != STATUS_CODE_NOT_ALLOWED)
                {
                    std::cout << "ElectronStartRefinementTool: cluster jam" << std::endl;
                    throw StatusCodeException(statusCodeCluster);
                }

                if (pParentPfo)
                {
                    const StatusCode statusCodePfo(PandoraContentApi::RemoveFromPfo(*this, pParentPfo, pParentCluster));
                    const unsigned int nHits(LArPfoHelper::GetNumberOfTwoDHits(pParentPfo));

                    if (nHits == 0)
                        std::cout << "ElectronStartRefinementTool: ISOBEL - PFO HAS ZERO HITS" << std::endl;

                    if (statusCodePfo != STATUS_CODE_SUCCESS)
                    {
                        std::cout << "ElectronStartRefinementTool: pfo jam" << std::endl;
                        throw;
                    }
                }

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pParentCluster));
            }
        }

        if (!PandoraContentApi::IsAvailable(*this, pPathwayHit))
        {
            std::cout << "CALO HIT IS NOT AVAILABLE!!" << std::endl;
            throw;
        }

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::AddToCluster(*this, showerClusters2D.front(), pPathwayHit));
    }

    // Now handle the shower merges
    for (const ParticleFlowObject *const pShowerToMerge : significantShowersToMerge)
    {
        m_deletedPfos.push_back(pShowerToMerge);
        this->MergeAndDeletePfos(pShowerPfo, pShowerToMerge);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HybridShowerStartRefinementAlgorithm::SetElectronMetadata(const CartesianVector &nuVertexPosition, const ParticleFlowObject *const pShowerPfo)
{
    object_creation::ParticleFlowObject::Metadata metadata;
    metadata.m_propertiesToAdd["ShowerVertexX"] = nuVertexPosition.GetX();
    metadata.m_propertiesToAdd["ShowerVertexY"] = nuVertexPosition.GetY();
    metadata.m_propertiesToAdd["ShowerVertexZ"] = nuVertexPosition.GetZ();
    metadata.m_propertiesToAdd["dEdX"] = 2.3;
    metadata.m_propertiesToAdd["ActiveHybridElectronAlg"] = 1.f;

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::AlterMetadata(*this, pShowerPfo, metadata));
}

//------------------------------------------------------------------------------------------------------------------------------------------

/*
void HybridShowerStartRefinementAlgorithm::RemoveConnectionPathway(const ParticleFlowObject *const pShowerPfo, const ProtoShower &protoShower)
{
    const HitType hitType(protoShower.m_connectionPathway.m_pathwayHitList.front()->GetHitType());

    ClusterList clusterList;
    LArPfoHelper::GetClusters(pShowerPfo, hitType, clusterList);

    std::string clusterListName(hitType == TPC_VIEW_U ? "ClustersU" : hitType == TPC_VIEW_V ? "ClustersV" : "ClustersW");
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));

    for (const CaloHit *const pCaloHit : protoShower.m_connectionPathway.m_pathwayHitList)
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::RemoveFromCluster(*this, clusterList.front(), pCaloHit));
}
*/
//------------------------------------------------------------------------------------------------------------------------------------------

void HybridShowerStartRefinementAlgorithm::RemoveConnectionPathway(const ParticleFlowObject *const pGammaPfo, const ElectronProtoShowerVector &protoShowerVectorU, 
    const ElectronProtoShowerVector &protoShowerVectorV, const ElectronProtoShowerVector &protoShowerVectorW)
{
    CaloHitList hitsToRemove;

    for (const ElectronProtoShower &protoShowerU : protoShowerVectorU)
    {
        if (!protoShowerU.m_connectionPathway.m_pathwayHitList.empty())
            hitsToRemove.insert(hitsToRemove.end(), protoShowerU.m_connectionPathway.m_pathwayHitList.begin(), protoShowerU.m_connectionPathway.m_pathwayHitList.end());
    }

    for (const ElectronProtoShower &protoShowerV : protoShowerVectorV)
    {
        if (!protoShowerV.m_connectionPathway.m_pathwayHitList.empty())
            hitsToRemove.insert(hitsToRemove.end(), protoShowerV.m_connectionPathway.m_pathwayHitList.begin(), protoShowerV.m_connectionPathway.m_pathwayHitList.end());
    }

    for (const ElectronProtoShower &protoShowerW : protoShowerVectorW)
    {
        if (!protoShowerW.m_connectionPathway.m_pathwayHitList.empty())
            hitsToRemove.insert(hitsToRemove.end(), protoShowerW.m_connectionPathway.m_pathwayHitList.begin(), protoShowerW.m_connectionPathway.m_pathwayHitList.end());
    }

    ClusterList twoDClusterList;
    LArPfoHelper::GetClusters(pGammaPfo, TPC_VIEW_U, twoDClusterList);
    LArPfoHelper::GetClusters(pGammaPfo, TPC_VIEW_V, twoDClusterList);
    LArPfoHelper::GetClusters(pGammaPfo, TPC_VIEW_W, twoDClusterList);

    MCParticleVector mcContaminantVector;
    MCParticleToHitListMap mcParticleToHitListMapU, mcParticleToHitListMapV, mcParticleToHitListMapW;

    for (const Cluster *const pCluster : twoDClusterList)
    {
        CaloHitList caloHitList;
        pCluster->GetOrderedCaloHitList().FillCaloHitList(caloHitList);

        // Add in isolated hits
        const CaloHitList &isolatedHitList(pCluster->GetIsolatedCaloHitList());

        for (const CaloHit *const pIsolated : isolatedHitList)
        {
            if (std::find(caloHitList.begin(), caloHitList.end(), pIsolated) == caloHitList.end())
                caloHitList.push_back(pIsolated);
        }

        const HitType hitType(caloHitList.front()->GetHitType());
        MCParticleToHitListMap &viewMCParticleToHitListMap(hitType == TPC_VIEW_U ? mcParticleToHitListMapU : hitType == TPC_VIEW_V ? mcParticleToHitListMapV : mcParticleToHitListMapW);
        const std::string &clusterListName(hitType == TPC_VIEW_U ? "ClustersU" : hitType == TPC_VIEW_V ? "ClustersV" : "ClustersW");

        CaloHitList removedIsolatedHits;

        // Remove Hits
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            if (std::find(hitsToRemove.begin(), hitsToRemove.end(), pCaloHit) == hitsToRemove.end())
                continue;

            if (std::find(removedIsolatedHits.begin(), removedIsolatedHits.end(), pCaloHit) != removedIsolatedHits.end())
                continue;

            const MCParticle *pMCParticle(nullptr);

            try
            {
                pMCParticle = MCParticleHelper::GetMainMCParticle(pCaloHit);
            }
            catch (const StatusCodeException &)
            {
                continue;
            }

            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, clusterListName));

            CaloHitList clusterNormalHitList;
            pCluster->GetOrderedCaloHitList().FillCaloHitList(clusterNormalHitList);
            const CaloHitList clusterIsolatedHitList(pCluster->GetIsolatedCaloHitList());

            const bool isIsolated(std::find(isolatedHitList.begin(), isolatedHitList.end(), pCaloHit) != isolatedHitList.end());

            if (!isIsolated && (clusterNormalHitList.size() == 1) && !(clusterIsolatedHitList.empty()))
            {
                for (const CaloHit * const pIsolatedHit : clusterIsolatedHitList)
                {
                    removedIsolatedHits.push_back(pIsolatedHit);
                    const StatusCode isolatedStatusCode(PandoraContentApi::RemoveIsolatedFromCluster(*this, pCluster, pIsolatedHit));

                    if (isolatedStatusCode != STATUS_CODE_SUCCESS)
                    {
                        std::cout << "ISOBEL CANNOT REMOVE ISOLATED HIT?" << std::endl;
                        throw;
                    }
                }
            }

            const StatusCode statusCodeCluster(isIsolated ? PandoraContentApi::RemoveIsolatedFromCluster(*this, pCluster, pCaloHit) :
                PandoraContentApi::RemoveFromCluster(*this, pCluster, pCaloHit));

            if (statusCodeCluster != STATUS_CODE_SUCCESS)
            {
                if (statusCodeCluster != STATUS_CODE_NOT_ALLOWED)
                {
                    std::cout << "HybridShowerStartRefinementAlgorithm: cluster jam" << std::endl;
                    throw StatusCodeException(statusCodeCluster);
                }

                const StatusCode statusCodePfo(PandoraContentApi::RemoveFromPfo(*this, pGammaPfo, pCluster));

                if (statusCodePfo != STATUS_CODE_SUCCESS)
                {
                    std::cout << "HybridShowerStartRefinementAlgorithm: pfo jam" << std::endl;
                    throw;
                }

                PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Delete(*this, pCluster));
            }

            if (!PandoraContentApi::IsAvailable(*this, pCaloHit))
            {
                std::cout << "HybridShowerStartRefinementAlgorithm: CALO HIT IS NOT AVAILABLE!!" << std::endl;
                throw;
            }
    
            viewMCParticleToHitListMap[pMCParticle].push_back(pCaloHit);

            if (std::find(mcContaminantVector.begin(), mcContaminantVector.end(), pMCParticle) == mcContaminantVector.end())
                mcContaminantVector.push_back(pMCParticle);
        }
    }

    // Do some particle creation

    std::sort(mcContaminantVector.begin(), mcContaminantVector.end(), LArMCParticleHelper::SortByMomentum);

    for (const MCParticle *const pMCParticle : mcContaminantVector)
    {
        PandoraContentApi::ParticleFlowObject::Parameters pfoParameters;
        pfoParameters.m_particleId = (this->IsShower(pMCParticle) ? 11 : 13);
        pfoParameters.m_charge = PdgTable::GetParticleCharge(pfoParameters.m_particleId.Get());
        pfoParameters.m_mass = PdgTable::GetParticleMass(pfoParameters.m_particleId.Get());
        pfoParameters.m_energy = pMCParticle->GetEnergy();
        pfoParameters.m_momentum = pMCParticle->GetMomentum();

        CaloHitList allHits(0);

        if (mcParticleToHitListMapU.find(pMCParticle) != mcParticleToHitListMapU.end())
            allHits.insert(allHits.end(), mcParticleToHitListMapU.at(pMCParticle).begin(), mcParticleToHitListMapU.at(pMCParticle).end());

        if (mcParticleToHitListMapV.find(pMCParticle) != mcParticleToHitListMapV.end())
            allHits.insert(allHits.end(), mcParticleToHitListMapV.at(pMCParticle).begin(), mcParticleToHitListMapV.at(pMCParticle).end());

        if (mcParticleToHitListMapW.find(pMCParticle) != mcParticleToHitListMapW.end())
            allHits.insert(allHits.end(), mcParticleToHitListMapW.at(pMCParticle).begin(), mcParticleToHitListMapW.at(pMCParticle).end());

        if (allHits.size() < 50)
            continue;

        if (mcParticleToHitListMapU.find(pMCParticle) != mcParticleToHitListMapU.end())
            pfoParameters.m_clusterList.push_back(this->CreateCluster(pMCParticle, mcParticleToHitListMapU.at(pMCParticle), TPC_VIEW_U));

        if (mcParticleToHitListMapV.find(pMCParticle) != mcParticleToHitListMapV.end())
            pfoParameters.m_clusterList.push_back(this->CreateCluster(pMCParticle, mcParticleToHitListMapV.at(pMCParticle), TPC_VIEW_V));

        if (mcParticleToHitListMapW.find(pMCParticle) != mcParticleToHitListMapW.end())
            pfoParameters.m_clusterList.push_back(this->CreateCluster(pMCParticle, mcParticleToHitListMapW.at(pMCParticle), TPC_VIEW_W));

        const PfoList *pTemporaryList(nullptr);
        std::string temporaryListName, currentListName;

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<ParticleFlowObject>(*this, currentListName));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
            PandoraContentApi::CreateTemporaryListAndSetCurrent<PfoList>(*this, pTemporaryList, temporaryListName));

        const ParticleFlowObject *pPfo(nullptr);

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ParticleFlowObject::Create(*this, pfoParameters, pPfo));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, 
            PandoraContentApi::SaveList<ParticleFlowObject>(*this, temporaryListName, this->IsShower(pMCParticle) ? "ShowerParticles3D" : "TrackParticles3D"));

        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<ParticleFlowObject>(*this, currentListName));
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool HybridShowerStartRefinementAlgorithm::IsShower(const MCParticle *const pMCParticle)
{
    const int pdg(pMCParticle->GetParticleId());
    return ((E_MINUS == std::abs(pdg)) || (PHOTON == std::abs(pdg)));
}

//------------------------------------------------------------------------------------------------------------------------------------------

const Cluster *HybridShowerStartRefinementAlgorithm::CreateCluster(const MCParticle *const pMCParticle, const CaloHitList &caloHitList, const HitType hitType)
{
    const ClusterList *pTemporaryList(nullptr);
    std::string temporaryListName, currentListName;

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetCurrentListName<Cluster>(*this, currentListName));

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=,
        PandoraContentApi::CreateTemporaryListAndSetCurrent<ClusterList>(*this, pTemporaryList, temporaryListName));

    const Cluster *pCluster(nullptr);
    PandoraContentApi::Cluster::Parameters parameters;
    parameters.m_caloHitList = caloHitList;
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::Create(*this, parameters, pCluster));

    PandoraContentApi::Cluster::Metadata metadata;
    metadata.m_particleId = this->IsShower(pMCParticle) ? 11 : 13;

    if (metadata.m_particleId.IsInitialized())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::Cluster::AlterMetadata(*this, pCluster, metadata));

    const std::string &clusterListName(hitType == TPC_VIEW_U ? "ClustersU" : hitType == TPC_VIEW_V ? "ClustersV" : "ClustersW");

    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList<Cluster>(*this, temporaryListName, clusterListName));
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::ReplaceCurrentList<Cluster>(*this, currentListName));

    return pCluster;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HybridShowerStartRefinementAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadVectorOfValues(xmlHandle, "PfoListNames", m_pfoListNames));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoVertexListName", m_neutrinoVertexListName));

    AlgorithmToolVector algorithmToolVector;
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmToolList(*this, xmlHandle, "ShowerStartRefinementTools", algorithmToolVector));

    for (AlgorithmToolVector::const_iterator iter = algorithmToolVector.begin(), iterEnd = algorithmToolVector.end(); iter != iterEnd; ++iter)
    {
        HybridShowerStartRefinementBaseTool *const pShowerStartRefinementTool(dynamic_cast<HybridShowerStartRefinementBaseTool *>(*iter));

        if (!pShowerStartRefinementTool)
            return STATUS_CODE_INVALID_PARAMETER;

        m_algorithmToolVector.push_back(pShowerStartRefinementTool);
    }

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, 
        XmlHelper::ReadValue(xmlHandle, "BinSize", m_binSize));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "ThresholdSignalGammaDisplacement", m_thresholdSignalGammaDisplacement));

    PfoMopUpBaseAlgorithm::ReadSettings(xmlHandle);

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content

    /*
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "ShowerDistribution", "ShowerDistribution.root", "UPDATE"));
    }
    catch (const StatusCodeException &)
    {
        std::cout << "THE LIMIT DOES NOT EXIST" << std::endl;
    }
    */