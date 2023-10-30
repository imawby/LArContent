/**
 *  @file   larpandoracontent/LArShowerRefinement/CPValidationAlgorithm.cc
 *
 *  @brief  Implementation of the electron initial region refinement algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingShowerStartFinderTool.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"
#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"
#include "larpandoracontent/LArHelpers/LArMvaHelper.h"
#include "larpandoracontent/LArHelpers/LArPfoHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingShowerFitResult.h"

#include "larpandoracontent/LArShowerRefinement/CPValidationAlgorithm.h"
#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"
#include "larpandoracontent/LArShowerRefinement/PeakDirectionFinderTool.h"
#include "larpandoracontent/LArShowerRefinement/ProtoShowerMatchingTool.h"
#include "larpandoracontent/LArShowerRefinement/ShowerSpineFinderTool.h"
#include "larpandoracontent/LArShowerRefinement/ShowerStartFinderTool.h"

#ifdef MONITORING
#include "PandoraMonitoringApi.h"
#endif

using namespace pandora;

namespace lar_content
{

CPValidationAlgorithm::CPValidationAlgorithm() :
    m_minShowerHits3D(50),
    m_showerSlidingFitWindow(1000),
    m_maxCoincidenceTransverseSeparation(5.f),
    m_minSpinePurity(0.7f),
    m_minElectronCompleteness(0.33f),
    m_minElectronPurity(0.5f),
    m_maxSeparationFromHit(3.f),
    m_maxProjectionSeparation(5.f),
    m_maxXSeparation(0.5f),
    m_cheatShowerStart(false),
    m_eventCounter(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

CPValidationAlgorithm::~CPValidationAlgorithm()
{
    try
    {
        PANDORA_MONITORING_API(SaveTree(this->GetPandora(), "CPValidationTree", "CPValidation.root", "UPDATE"));
    }
    catch (const StatusCodeException &)
    {
        std::cout << "CPValidationAlgorithm: Unable to write tree CPValidationTree to file CPValidation.root" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CPValidationAlgorithm::Run()
{
    ++m_eventCounter;
    //////////////////////////////////////////
    PANDORA_MONITORING_API(SetEveDisplayParameters(this->GetPandora(), true, DETECTOR_VIEW_XZ, -1.f, -1.f, 1.f));
    //////////////////////////////////////////

    PfoVector showerPfoVector;
    this->FillShowerPfoVector(showerPfoVector);

    //////////////////////////////////////////   
    //std::cout << "showerPfoVector.size(): " << showerPfoVector.size() << std::endl;
    //////////////////////////////////////////   

    if (showerPfoVector.empty())
        return STATUS_CODE_SUCCESS;

    for (const ParticleFlowObject *const pShowerPfo : showerPfoVector)
    {
        m_foundCPU = 0;
        m_foundCPV = 0;
        m_foundCPW = 0;
        m_matchIn3D = 0;
        m_createdShowerStarts3D = 0;

        this->RefineShower(pShowerPfo);

        //////////////////////////////////////////  
        /*
        std::cout << "FoundCPU: " << (m_foundCPU ? "yes" : "no") << std::endl;
        std::cout << "FoundCPV: " << (m_foundCPV ? "yes" : "no") << std::endl;
        std::cout << "FoundCPW: " << (m_foundCPW ? "yes" : "no") << std::endl;
        std::cout << "MatchIn3D: " << (m_matchIn3D ? "yes" : "no") << std::endl;
        std::cout << "CreatedShowerStarts3D: " << (m_createdShowerStarts3D ? "yes" : "no") << std::endl;
        PandoraMonitoringApi::ViewEvent(this->GetPandora()); 
        */
        //////////////////////////////////////////

        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CPValidationTree", "EventCount", m_eventCounter));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CPValidationTree", "FoundCPU", m_foundCPU));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CPValidationTree", "FoundCPV", m_foundCPV));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CPValidationTree", "FoundCPW", m_foundCPW));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CPValidationTree", "MatchIn3D", m_matchIn3D));
        PANDORA_MONITORING_API(SetTreeVariable(this->GetPandora(), "CPValidationTree", "CreatedShowerStarts3D", m_createdShowerStarts3D));

        PANDORA_MONITORING_API(FillTree(this->GetPandora(), "CPValidationTree"));
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CPValidationAlgorithm::FillShowerPfoVector(PfoVector &showerPfoVector) const
{
    // Get true electron hit map
    HitOwnershipMap electronHitMap;
    this->FillElectronHitMap(electronHitMap);

    const PfoList *pPfoList(nullptr);

    if (PandoraContentApi::GetList(*this, m_showerPfoListName, pPfoList) != STATUS_CODE_SUCCESS)
        return;

    if (!pPfoList || pPfoList->empty())
    {
        std::cout << "CPValidationAlgorithm: unable to find shower pfo list " << m_showerPfoListName << std::endl;
        return;
    }

    for (const ParticleFlowObject *const pShowerPfo : *pPfoList)
    {
        // I only care about the CPs of true primary electrons
        if (!this->IsElectron(pShowerPfo, electronHitMap))
            continue;

        // Only consider significant showers
        CaloHitList caloHits3D;
        LArPfoHelper::GetCaloHits(pShowerPfo, TPC_3D, caloHits3D);

        if (caloHits3D.size() < m_minShowerHits3D)
            continue;

        showerPfoVector.push_back(pShowerPfo);
    }

    std::sort(showerPfoVector.begin(), showerPfoVector.end(), LArPfoHelper::SortByNHits);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CPValidationAlgorithm::RefineShower(const ParticleFlowObject *const pShowerPfo)
{
    CartesianVector nuVertex3D(0.f, 0.f, 0.f);

    if (this->GetNeutrinoVertex(nuVertex3D) != STATUS_CODE_SUCCESS)
        return;

    // Create the 2D connetion pathways
    ProtoShowerVector protoShowerVectorU, protoShowerVectorV, protoShowerVectorW;

    this->BuildViewProtoShowers(pShowerPfo, nuVertex3D, TPC_VIEW_U, protoShowerVectorU);

    if (!protoShowerVectorU.empty())
        m_foundCPU = 1;

    this->BuildViewProtoShowers(pShowerPfo, nuVertex3D, TPC_VIEW_V, protoShowerVectorV);

    if (!protoShowerVectorV.empty())
        m_foundCPV = 1;

    this->BuildViewProtoShowers(pShowerPfo, nuVertex3D, TPC_VIEW_W, protoShowerVectorW);

    if (!protoShowerVectorW.empty())
        m_foundCPW = 1;

    /////////////////////////////////////
    /*
    for (const ProtoShower &protoShower : protoShowerVectorU)
    {
        CaloHitList showerCaloHitsU;
        LArPfoHelper::GetCaloHits(pShowerPfo, TPC_VIEW_U, showerCaloHitsU);
        PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &showerCaloHitsU, "showerCaloHitsU", RED);
        PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &protoShower.GetSpineHitList(), "spineHitListU", VIOLET);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &protoShower.GetShowerCore().GetStartPosition(), "showerStartU", BLUE, 2);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }

    for (const ProtoShower &protoShower : protoShowerVectorV)
    {
        CaloHitList showerCaloHitsV;
        LArPfoHelper::GetCaloHits(pShowerPfo, TPC_VIEW_V, showerCaloHitsV);
        PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &showerCaloHitsV, "showerCaloHitsV", RED);
        PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &protoShower.GetSpineHitList(), "spineHitListV", VIOLET);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &protoShower.GetShowerCore().GetStartPosition(), "showerStartV", BLUE, 2);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }

    for (const ProtoShower &protoShower : protoShowerVectorW)
    {
        CaloHitList showerCaloHitsW;
        LArPfoHelper::GetCaloHits(pShowerPfo, TPC_VIEW_W, showerCaloHitsW);
        PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &showerCaloHitsW, "showerCaloHitsW", RED);
        PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &protoShower.GetSpineHitList(), "spineHitListW", VIOLET);
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &protoShower.GetShowerCore().GetStartPosition(), "showerStartW", BLUE, 2);
        PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }
    */
    /////////////////////////////////////

    // 2D->3D connection pathway matching
    ProtoShowerMatchVector protoShowerMatchVector;
    m_pProtoShowerMatchingTool->Run(protoShowerVectorU, protoShowerVectorV, protoShowerVectorW, protoShowerMatchVector);

    if (!protoShowerMatchVector.empty())
        m_matchIn3D = 1;

    for (ProtoShowerMatch &protoShowerMatch : protoShowerMatchVector)
    {
        if (!m_createdShowerStarts3D)
        {
            // Determine the 3D shower vertex
            CartesianPointVector showerStarts3D;

            if (LArConnectionPathwayHelper::FindShowerStarts3D(this, pShowerPfo, protoShowerMatch, nuVertex3D, m_maxSeparationFromHit,
                m_maxProjectionSeparation, m_maxXSeparation, showerStarts3D))
            {
                m_createdShowerStarts3D = 1;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CPValidationAlgorithm::GetNeutrinoVertex(CartesianVector &nuVertex3D) const
{
    const VertexList *pNuVertexList(nullptr);
    const StatusCode statusCode(PandoraContentApi::GetList(*this, m_neutrinoVertexListName, pNuVertexList));

    if (statusCode != STATUS_CODE_SUCCESS)
        return statusCode;

    if (!pNuVertexList || (pNuVertexList->size() != 1))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "CPValidationAlgorithm: unable to find vertex list " << m_neutrinoVertexListName
                      << " if it does exist, it may have more than one nu vertex" << std::endl;

        return STATUS_CODE_NOT_INITIALIZED;
    }

    nuVertex3D = pNuVertexList->front()->GetPosition();

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CPValidationAlgorithm::BuildViewProtoShowers(const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertex3D, HitType hitType, ProtoShowerVector &protoShowerVector) const
{
    const CaloHitList *pViewHitList(nullptr);

    if (this->GetHitListOfType(hitType, pViewHitList) != STATUS_CODE_SUCCESS)
        return;

    CartesianVector showerVertexPosition(0.f, 0.f, 0.f);
    try
    {
        showerVertexPosition = this->GetShowerVertex(pShowerPfo, hitType, nuVertex3D);
    }
    catch (...)
    {
        return;
    }

    // Determine directions of pathways out of neutrino vertex
    CartesianPointVector peakDirectionVector;
    if (m_pShowerPeakDirectionFinderTool->Run(pShowerPfo, nuVertex3D, pViewHitList, hitType, peakDirectionVector) != STATUS_CODE_SUCCESS)
        return;

    // Investigate each direction
    CaloHitList unavailableHitList;
    for (CartesianVector &peakDirection : peakDirectionVector)
    {
        // Create the ProtoShower object
        const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));

        //const CartesianVector end(nuVertex2D + (peakDirection * 20.f));
        //PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &nuVertex2D, &end, "FU", BLACK, 2, 2);

        // Collect the hits associated with the pathway (the shower spine)
        CaloHitList showerSpineHitList;
        if (m_pShowerSpineFinderTool->Run(nuVertex3D, pViewHitList, hitType, peakDirection, unavailableHitList, showerSpineHitList) != STATUS_CODE_SUCCESS)
            continue;

        //std::cout << "A" << std::endl;

        this->RefineShowerVertex(pShowerPfo, hitType, nuVertex3D, peakDirection, showerVertexPosition);

        //std::cout << "B" << std::endl;

        // Quality check: If the spine passes the shower vertex, does it live inside the shower?
        if (!this->IsSpineCoincident(pShowerPfo, nuVertex3D, hitType, showerVertexPosition, showerSpineHitList))
            continue;

        //std::cout << "C" << std::endl;

        // Find the 2D shower start position
        CartesianVector showerStartPosition(0.f, 0.f, 0.f);

        if (m_cheatShowerStart)
        {
            if (!this->GetShowerStart(pShowerPfo, hitType, showerSpineHitList, showerStartPosition))
                continue;
        }
        else
        {
            CartesianVector showerStartDirection(0.f, 0.f, 0.f);

            if (m_pShowerStartFinderTool->Run(pShowerPfo, peakDirection, hitType, showerSpineHitList, showerStartPosition, showerStartDirection) != STATUS_CODE_SUCCESS)
                continue;
        }

        //std::cout << "D" << std::endl;

        ProtoShower protoShower(ShowerCore(showerStartPosition), ConnectionPathway(nuVertex2D, peakDirection),
            showerSpineHitList, CaloHitList(), CartesianPointVector(), CaloHitList());

        // Now determine the hits to be added to the shower from the found spine
        CaloHitList viewShowerHitList;
        LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewShowerHitList);

        const float showerVertexL(std::max(
            peakDirection.GetDotProduct(showerStartPosition - nuVertex2D), peakDirection.GetDotProduct(showerVertexPosition - nuVertex2D)));

        for (const CaloHit *const pCaloHit : showerSpineHitList)
        {
            if (std::find(viewShowerHitList.begin(), viewShowerHitList.end(), pCaloHit) != viewShowerHitList.end())
                continue;

            const CartesianVector &hitPosition(pCaloHit->GetPositionVector());
            const float showerL(peakDirection.GetDotProduct(hitPosition - nuVertex2D));

            if ((showerL > 0.f) && (showerL < showerVertexL))
                protoShower.AddHitToAdd(pCaloHit);
        }

        protoShowerVector.push_back(protoShower);
        unavailableHitList.insert(unavailableHitList.begin(), showerSpineHitList.begin(), showerSpineHitList.end());
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CPValidationAlgorithm::GetHitListOfType(const HitType hitType, const CaloHitList *&pCaloHitList) const
{
    const std::string typeHitListName(hitType == TPC_VIEW_U ? m_caloHitListNameU : hitType == TPC_VIEW_V ? m_caloHitListNameV : m_caloHitListNameW);

    PANDORA_THROW_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_INITIALIZED, !=, PandoraContentApi::GetList(*this, typeHitListName, pCaloHitList));

    if (!pCaloHitList || pCaloHitList->empty())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "ShowerStartRefinementBaseTool: unable to find calo hit list " << typeHitListName << std::endl;

        return STATUS_CODE_NOT_INITIALIZED;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

CartesianVector CPValidationAlgorithm::GetShowerVertex(
    const ParticleFlowObject *const pShowerPfo, const HitType hitType, const CartesianVector &nuVertex3D) const
{
    ClusterList viewCusterList;
    LArPfoHelper::GetClusters(pShowerPfo, hitType, viewCusterList);

    if (viewCusterList.empty())
        throw StatusCodeException(STATUS_CODE_FAILURE);

    const TwoDSlidingShowerFitResult twoDShowerSlidingFit(
        viewCusterList.front(), m_showerSlidingFitWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), hitType));
    const CartesianVector &minLayerPosition(twoDShowerSlidingFit.GetShowerFitResult().GetGlobalMinLayerPosition());
    const CartesianVector &maxLayerPosition(twoDShowerSlidingFit.GetShowerFitResult().GetGlobalMaxLayerPosition());
    const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));

    if ((nuVertex2D.GetZ() > minLayerPosition.GetZ()) && (nuVertex2D.GetZ() < maxLayerPosition.GetZ()))
        return nuVertex2D;

    const float minSeparation((nuVertex2D - minLayerPosition).GetMagnitudeSquared());
    const float maxSeparation((nuVertex2D - maxLayerPosition).GetMagnitudeSquared());
    CartesianVector showerVertexPosition(minSeparation < maxSeparation ? minLayerPosition : maxLayerPosition);

    return (minSeparation < maxSeparation ? minLayerPosition : maxLayerPosition);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CPValidationAlgorithm::RefineShowerVertex(const ParticleFlowObject *const pShowerPfo, const HitType hitType,
    const CartesianVector &nuVertex3D, const CartesianVector &peakDirection, CartesianVector &showerVertexPosition) const
{
    const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));

    if (!this->IsShowerConnected(showerVertexPosition, nuVertex2D, peakDirection))
    {
        CaloHitList viewShowerHitList;
        LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewShowerHitList);

        float minL(std::numeric_limits<float>::max());

        for (const CaloHit *const pCaloHit : viewShowerHitList)
        {
            const CartesianVector displacement(pCaloHit->GetPositionVector() - nuVertex2D);
            const float longitudinalSeparation(peakDirection.GetDotProduct(displacement));

            if ((longitudinalSeparation < (-1.f)) || (longitudinalSeparation > minL))
                continue;

            const float transverseSeparation((peakDirection.GetCrossProduct(displacement)).GetMagnitude());

            if (transverseSeparation < m_maxCoincidenceTransverseSeparation)
            {
                showerVertexPosition = pCaloHit->GetPositionVector();
                minL = longitudinalSeparation;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CPValidationAlgorithm::IsShowerConnected(
    const CartesianVector &showerVertexPosition, const CartesianVector &nuVertex2D, const CartesianVector &peakDirection) const
{
    CartesianVector displacement(showerVertexPosition - nuVertex2D);

    const float longitudinalSeparation(peakDirection.GetDotProduct(displacement));

    if (longitudinalSeparation < (-1.f))
        return false;

    const float transverseSeparation((peakDirection.GetCrossProduct(displacement)).GetMagnitude());

    if (transverseSeparation > m_maxCoincidenceTransverseSeparation)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CPValidationAlgorithm::IsSpineCoincident(const ParticleFlowObject *const pShowerPfo,
    const CartesianVector &nuVertex3D, const HitType hitType, const CartesianVector &showerVertex, const CaloHitList &showerSpineHitList) const
{
    const CartesianVector nuVertex2D(LArGeometryHelper::ProjectPosition(this->GetPandora(), nuVertex3D, hitType));

    CaloHitList viewShowerHitList;
    LArPfoHelper::GetCaloHits(pShowerPfo, hitType, viewShowerHitList);

    CaloHitList postShowerVertexSpineHits;
    const float showerDistanceFromNuSquared((showerVertex - nuVertex2D).GetMagnitudeSquared());

    for (const CaloHit *const pSpineHit : showerSpineHitList)
    {
        const CartesianVector &hitPosition(pSpineHit->GetPositionVector());
        const float separationSquared((hitPosition - nuVertex2D).GetMagnitudeSquared());

        if (separationSquared > showerDistanceFromNuSquared)
            postShowerVertexSpineHits.push_back(pSpineHit);
    }

    if (postShowerVertexSpineHits.size() == 0)
        return true;

    // Check whether shower spine is pure
    const CaloHitList sharedHitList(LArMCParticleHelper::GetSharedHits(postShowerVertexSpineHits, viewShowerHitList));
    const float spinePurity(static_cast<float>(sharedHitList.size()) / static_cast<float>(postShowerVertexSpineHits.size()));

    //std::cout << "spinePurity: " << spinePurity << std::endl;

    if (spinePurity < m_minSpinePurity)
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CPValidationAlgorithm::GetShowerStart(const ParticleFlowObject *const pPfo, const HitType hitType, 
    const CaloHitList &showerSpineHitList, CartesianVector &showerStart) const
{
    // Combine HitLists
    CaloHitList combinedHitList;
    LArPfoHelper::GetCaloHits(pPfo, hitType, combinedHitList);

    for (const CaloHit *const pCaloHit : showerSpineHitList)
    {
        if (std::find(combinedHitList.begin(), combinedHitList.end(), pCaloHit) == combinedHitList.end())
            combinedHitList.push_back(pCaloHit);
    }

    // Get MCParticle List
    const MCParticleList *pMCParticleList(nullptr);

    if (PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList) != STATUS_CODE_SUCCESS)
        return false;

    if (!pMCParticleList || pMCParticleList->empty())
    {
        std::cout << "MLShowerInitialRegionRefinementAlgorithm: unable to find mc particle list " << m_mcParticleListName << std::endl;
        return false;
    }

    return (m_pCheatingShowerStartFinderTool->Run(pPfo, pMCParticleList, combinedHitList, hitType, showerStart) == STATUS_CODE_SUCCESS);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

void CPValidationAlgorithm::FillElectronHitMap(HitOwnershipMap &electronHitMap) const
{
    electronHitMap.clear();

    for (HitType hitType : {TPC_VIEW_U, TPC_VIEW_V, TPC_VIEW_W})
    {
        const CaloHitList *pViewHitList(nullptr);

        if (this->GetHitListOfType(hitType, pViewHitList) != STATUS_CODE_SUCCESS)
            continue;

        for (const CaloHit *const pCaloHit : *pViewHitList)
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
                const bool isLeadingElectron((std::abs(pMCParticle->GetParticleId()) == 11) &&
                                             (LArMCParticleHelper::GetPrimaryMCParticle(pMCParticle) == pMCParticle));

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
                electronHitMap[highestElectronContributor].push_back(pCaloHit);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool CPValidationAlgorithm::IsElectron(const ParticleFlowObject *const pPfo, const HitOwnershipMap &electronHitMap) const
{
    MCParticleVector mcElectronVector;

    for (auto &entry : electronHitMap)
        mcElectronVector.push_back(entry.first);

    CaloHitList pfoHitList;
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_U, pfoHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_V, pfoHitList);
    LArPfoHelper::GetCaloHits(pPfo, TPC_VIEW_W, pfoHitList);

    for (const MCParticle *const pMCElectron : mcElectronVector)
    {
        const CaloHitList &mcElectronHitList(electronHitMap.at(pMCElectron));
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

StatusCode CPValidationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "ShowerPfoListName", m_showerPfoListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "NeutrinoVertexListName", m_neutrinoVertexListName));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameU", m_caloHitListNameU));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameV", m_caloHitListNameV));
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListNameW", m_caloHitListNameW));

    AlgorithmTool *pAlgorithmTool1(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "ShowerPeakDirectionFinder", pAlgorithmTool1));
    m_pShowerPeakDirectionFinderTool = dynamic_cast<PeakDirectionFinderTool *>(pAlgorithmTool1);

    if (!m_pShowerPeakDirectionFinderTool)
        return STATUS_CODE_INVALID_PARAMETER;

    AlgorithmTool *pAlgorithmTool3(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "ShowerSpineFinder", pAlgorithmTool3));
    m_pShowerSpineFinderTool = dynamic_cast<ShowerSpineFinderTool *>(pAlgorithmTool3);

    if (!m_pShowerSpineFinderTool)
        return STATUS_CODE_INVALID_PARAMETER;

    AlgorithmTool *pAlgorithmTool5(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "ShowerStartFinder", pAlgorithmTool5));
    m_pShowerStartFinderTool = dynamic_cast<ShowerStartFinderTool *>(pAlgorithmTool5);

    if (!m_pShowerSpineFinderTool)
        return STATUS_CODE_INVALID_PARAMETER;

    AlgorithmTool *pAlgorithmTool6(nullptr);
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ProcessAlgorithmTool(*this, xmlHandle, "ProtoShowerMatching", pAlgorithmTool6));
    m_pProtoShowerMatchingTool = dynamic_cast<ProtoShowerMatchingTool *>(pAlgorithmTool6);

    if (!m_pProtoShowerMatchingTool)
        return STATUS_CODE_INVALID_PARAMETER;

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinShowerHits3D", m_minShowerHits3D));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "ShowerSlidingFitWindow", m_showerSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxCoincidenceTransverseSeparation", m_maxCoincidenceTransverseSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinSpinePurity", m_minSpinePurity));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MinElectronCompleteness", m_minElectronCompleteness));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSeparationFromHit", m_maxSeparationFromHit));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=,
        XmlHelper::ReadValue(xmlHandle, "MaxProjectionSeparation", m_maxProjectionSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxXSeparation", m_maxXSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinElectronPurity", m_minElectronPurity));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "CheatShowerStart", m_cheatShowerStart));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
