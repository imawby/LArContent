/**
 *  @file   larpandoradlcontent/LArThreeDReco/LArEventBuilding/DLSecondaryVertexSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the DL neutrino hierarchy algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoradlcontent/LArTwoDReco/DLSecondaryVertexSplittingAlgorithm.h"

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DLSecondaryVertexSplittingAlgorithm::DLSecondaryVertexSplittingAlgorithm() :
    m_visualise(false),
    m_clusterListNameU("ClustersU"),
    m_clusterListNameV("ClustersV"),
    m_clusterListNameW("ClustersW"),
    m_secondaryVertexListName("SecondaryVertices3D"),
    m_proximityForMatch(2.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLSecondaryVertexSplittingAlgorithm::Run()
{
    // Visualise
    if (m_visualise)
        PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, 1.f, 1.f));

    // Get secondary vertices
    const VertexList *pSecVtxList(nullptr);

    if (PandoraContentApi::GetList(*this, m_secondaryVertexListName, pSecVtxList) != STATUS_CODE_SUCCESS)
        return STATUS_CODE_SUCCESS;

    for (const HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        std::string clusterListName(hitType == TPC_VIEW_U ? m_clusterListNameU : hitType == TPC_VIEW_V ? m_clusterListNameV : m_clusterListNameW);
        std::string reclusterListName(hitType == TPC_VIEW_U ? "ReclusterU" : hitType == TPC_VIEW_V ? "ReclusterV" : "ReclusterW");

        // Get clusters
        const ClusterList *pClusterList(nullptr);

        if (PandoraContentApi::GetList(*this, clusterListName, pClusterList) != STATUS_CODE_SUCCESS)
            continue;

        // Create reclustering list
        ClusterList toRecluster;

        // Loop over clusters, searching for 'kink points'
        for (const Cluster *const pCluster : *pClusterList)
        {
            if (pCluster->GetNCaloHits() < 5)
                continue;

            try
            {
                const float slidingFitPitch(LArGeometryHelper::GetWirePitch(this->GetPandora(), LArClusterHelper::GetClusterHitType(pCluster)));
                const TwoDSlidingFitResult slidingFitResult(pCluster, 20, slidingFitPitch);

                // Get number of kink points along track
                const int nKinkPoints(this->GetNKinkPoints(slidingFitResult, pSecVtxList));

                if (nKinkPoints >= 1)
                    toRecluster.emplace_back(pCluster);
            }
            catch (...) { continue; }
        }

        if (!toRecluster.empty())
            PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, clusterListName, reclusterListName, toRecluster));
    }


    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

int DLSecondaryVertexSplittingAlgorithm::GetNKinkPoints(const TwoDSlidingFitResult &slidingFitResult, const VertexList *const pSecVtxList)
{
    int nKinkPoints(0);

    const Cluster *const pCluster(slidingFitResult.GetCluster());
    const HitType hitType(pCluster->GetOrderedCaloHitList().begin()->second->front()->GetHitType());
    const CartesianVector minPos(slidingFitResult.GetGlobalMinLayerPosition());
    const CartesianVector maxPos(slidingFitResult.GetGlobalMaxLayerPosition());

    for (const Vertex *const pSecVtx : *pSecVtxList)
    {
        const CartesianVector twoDPos(LArGeometryHelper::ProjectPosition(this->GetPandora(), pSecVtx->GetPosition(), hitType));

        if ((minPos - twoDPos).GetMagnitudeSquared() < (m_proximityForMatch * m_proximityForMatch))
            continue;

        if ((maxPos - twoDPos).GetMagnitudeSquared() < (m_proximityForMatch * m_proximityForMatch))
            continue;

        const float sep(LArClusterHelper::GetClosestDistance(twoDPos, pCluster));

        if (sep < m_proximityForMatch)
            ++nKinkPoints;

        //std::cout << "sep: " << sep << std::endl;
    }

    return nKinkPoints;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DLSecondaryVertexSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "Visualise", m_visualise));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameU", m_clusterListNameU));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameV", m_clusterListNameV));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ClusterListNameW", m_clusterListNameW));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SecondaryVertexListName", m_secondaryVertexListName));
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ProximityForMatch", m_proximityForMatch));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithm(*this, xmlHandle, "ClusterRebuilding", m_reclusteringAlgorithmName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_dl_content
