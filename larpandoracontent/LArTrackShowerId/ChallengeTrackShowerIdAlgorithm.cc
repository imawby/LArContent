/**
 *  @file   larpandoracontent/LArCheating/ChallengeTrackShowerIdAlgorithm.cc
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

#include "larpandoracontent/LArTrackShowerId/ChallengeTrackShowerIdAlgorithm.h"

using namespace pandora;

namespace lar_content
{

ChallengeTrackShowerIdAlgorithm::ChallengeTrackShowerIdAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ChallengeTrackShowerIdAlgorithm::Run()
{

    // Get pfos
    const PfoList *pTrackPfos(nullptr), *pShowerPfos(nullptr);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "TrackParticles3D", pTrackPfos));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "ShowerParticles3D", pShowerPfos));

    for (const PfoList *const pPfoList : {pTrackPfos, pShowerPfos})
    {
        if (!pPfoList || pPfoList->empty())
            continue;

        for (const ParticleFlowObject *const pPfo : *pPfoList)
        {
            ClusterList clustersU, clustersV, clustersW;
            LArPfoHelper::GetClusters(pPfo, TPC_VIEW_U, clustersU);
            LArPfoHelper::GetClusters(pPfo, TPC_VIEW_V, clustersV);
            LArPfoHelper::GetClusters(pPfo, TPC_VIEW_W, clustersW);

            int nU(clustersU.size()), nV(clustersV.size()), nW(clustersW.size());

            // Only look at two-view clusters
            if ((nU + nV + nW) != 2)
                continue;

            // Work out which views the clusters are in
            HitType hitType1(TPC_3D), hitType2(TPC_3D);
            for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
            {
                HitType &hitTypeToSet(hitType1 == TPC_3D ? hitType1 : hitType2);
                
                if (hitType == TPC_VIEW_U && nU == 1)
                    hitTypeToSet = TPC_VIEW_U;
                else if (hitType == TPC_VIEW_V && nV == 1)
                    hitTypeToSet = TPC_VIEW_V;
                else if (hitType == TPC_VIEW_W && nW == 1)
                    hitTypeToSet = TPC_VIEW_W;
            }
            
            const Cluster *const pCluster1(hitType1 == TPC_VIEW_U ? clustersU.front() : (hitType1 == TPC_VIEW_V ? clustersV.front() : clustersW.front()));
            const Cluster *const pCluster2(hitType2 == TPC_VIEW_U ? clustersU.front() : (hitType2 == TPC_VIEW_V ? clustersV.front() : clustersW.front()));

            // Now fit them both...
            try
            {
                const float slidingFitPitch1(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pCluster1)));
                const TwoDSlidingFitResult slidingFitResult1(pCluster1, 20, slidingFitPitch1);

                const float slidingFitPitch2(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pCluster2)));
                const TwoDSlidingFitResult slidingFitResult2(pCluster2, 20, slidingFitPitch2);

                float xMin1(0.f), xMax1(0.f), xMin2(0.f), xMax2(0.f);
                pCluster1->GetClusterSpanX(xMin1, xMax1);
                pCluster2->GetClusterSpanX(xMin2, xMax2);
                float xMin(std::max(xMin1, xMin2)), xMax(std::min(xMax1, xMax2));

                const int nSamplingPoints(std::floor(xMax - xMin) / 1.0);

                for (int i = 0; i < nSamplingPoints; ++i)
                {
                    const float xSample(xMin + i);
                    CartesianVector pos1(0.f, 0.f, 0.f), pos2(0.f, 0.f, 0.f);
                         
                    if (STATUS_CODE_SUCCESS != slidingFitResult1.GetGlobalFitPositionAtX(xSample, pos1))
                        continue;

                    if (STATUS_CODE_SUCCESS != slidingFitResult2.GetGlobalFitPositionAtX(xSample, pos2))
                        continue;
                         
                    const float z(LArGeometryHelper::MergeTwoPositions(this->GetPandora(), hitType1, hitType2, pos1.GetZ(), pos2.GetZ()));
                    CartesianVector proj(xSample, 0.f, z);
                    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &proj, "proj", BLUE, 2);
                }
            }
            catch (...)
            {
                continue;
            }            
        }        
        //PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &vis, "pTrackPfo", RED);
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }


    
    // Get u-view clusters
    const ClusterList *pClusterListU(nullptr);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "ClustersU", pClusterListU));

    ClusterList visU;
    for (const Cluster *const pCluster : *pClusterListU)
    {
        if (pCluster->IsAvailable())
            visU.push_back(pCluster);
    }

    // Get v-view clusters
    const ClusterList *pClusterListV(nullptr);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "ClustersV", pClusterListV));

    ClusterList visV;
    for (const Cluster *const pCluster : *pClusterListV)
    {
        if (pCluster->IsAvailable())
            visV.push_back(pCluster);
    }

    // Get w-view clusters
    const ClusterList *pClusterListW(nullptr);
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "ClustersW", pClusterListW));

    ClusterList visW;
    for (const Cluster *const pCluster : *pClusterListW)
    {
        if (pCluster->IsAvailable())
            visW.push_back(pCluster);
    }

    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &visU, "U-view clusters", RED);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &visV, "V-view clusters", GREEN);
    PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &visW, "W-view clusters", BLUE);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());

    
    
    // // Get shower pfos
    // const PfoList *pShowerPfoList(nullptr);
    // PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, "ShowerParticles3D", pShowerPfoList));

    // if (!pShowerPfoList || pShowerPfoList->empty())
    // {
    //     //if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
    //     //std::cout << "ChallengeTrackShowerIdAlgorithm: unable to find pfo list " << pfoListName << std::endl;

    //     return STATUS_CODE_SUCCESS;
    // }

    // for (const ParticleFlowObject *const pShowerPfo : *pShowerPfoList)
    // {
    //     const PropertiesMap &metadata(pShowerPfo->GetPropertiesMap());

    //     if (metadata.find("muon_pid_score") == metadata.end())
    //         continue;

    //     PfoList vis({pShowerPfo});
        
    //     std::cout << "Muon score: " << metadata.at("muon_pid_score") << std::endl;
    //     std::cout << "Proton score: " << metadata.at("proton_pid_score") << std::endl;
    //     std::cout << "Pion score: " << metadata.at("pion_pid_score") << std::endl;
    //     std::cout << "kaon score: " << metadata.at("kaon_pid_score") << std::endl;
    //     std::cout << "-------------------------------" << std::endl;

    //     PandoraMonitoringApi::VisualizeParticleFlowObjects(this->GetPandora(), &vis, "pShowerPfo", BLUE);
    //     PandoraMonitoringApi::ViewEvent(this->GetPandora());
    // }
    
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ChallengeTrackShowerIdAlgorithm::ReadSettings([[maybe_unused]] const TiXmlHandle xmlHandle)
{

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
