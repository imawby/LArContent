/**
 *  @file   larpandoracontent/LArShowerRefinement/ProtoShowerMatchingTool.cc
 *
 *  @brief  Implementation of the ProtoShower matching tool class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Pandora/AlgorithmTool.h"

#include "larpandoracontent/LArHelpers/LArConnectionPathwayHelper.h"
#include "larpandoracontent/LArHelpers/LArGeometryHelper.h"

#include "larpandoracontent/LArObjects/LArTwoDSlidingFitResult.h"

#include "larpandoracontent/LArShowerRefinement/LArProtoShower.h"
#include "larpandoracontent/LArShowerRefinement/ProtoShowerMatchingTool.h"

using namespace pandora;

namespace lar_content
{

ProtoShowerMatchingTool::ProtoShowerMatchingTool() :
    m_spineSlidingFitWindow(20),
    m_maxXSeparation(5.f),
    m_maxSeparation(5.f),
    m_maxAngularDeviation(5.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProtoShowerMatchingTool::Run(const ProtoShowerVector &protoShowerVectorU, const ProtoShowerVector &protoShowerVectorV,
    const ProtoShowerVector &protoShowerVectorW, ProtoShowerMatchVector &protoShowerMatchVector)
{
    IntVector usedProtoShowersU, usedProtoShowersV, usedProtoShowersW;

    bool added(false);

    for (unsigned int uIndex = 0; uIndex < protoShowerVectorU.size(); ++uIndex)
    {
        added = false;

        if (std::find(usedProtoShowersU.begin(), usedProtoShowersU.end(), uIndex) != usedProtoShowersU.end())
            continue;

        const ProtoShower &protoShowerU(protoShowerVectorU.at(uIndex));

        for (unsigned int vIndex = 0; vIndex < protoShowerVectorV.size(); ++vIndex)
        {
            if (std::find(usedProtoShowersV.begin(), usedProtoShowersV.end(), vIndex) != usedProtoShowersV.end())
                continue;

            const ProtoShower &protoShowerV(protoShowerVectorV.at(vIndex));

            for (unsigned int wIndex = 0; wIndex < protoShowerVectorW.size(); ++wIndex)
            {
                if (std::find(usedProtoShowersW.begin(), usedProtoShowersW.end(), wIndex) != usedProtoShowersW.end())
                    continue;

                const ProtoShower &protoShowerW(protoShowerVectorW.at(wIndex));
                Consistency consistency(Consistency::POSITION);

                if (!this->ArePathwaysConsistent(protoShowerU, protoShowerV, protoShowerW, consistency))
                    continue;

                usedProtoShowersU.push_back(uIndex);
                usedProtoShowersV.push_back(vIndex);
                usedProtoShowersW.push_back(wIndex);

                protoShowerMatchVector.push_back(ProtoShowerMatch(protoShowerU, protoShowerV, protoShowerW, consistency));

                added = true;
                break;
            }

            if (added)
                break;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ProtoShowerMatchingTool::ArePathwaysConsistent(
    const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW, Consistency &consistency) const
{
    //////////////////////////
    /*
    std::cout << "Are Pathways Consistent?" << std::endl;

    PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &protoShowerU.GetSpineHitList(), "spineHitListU", VIOLET);
    PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &protoShowerV.GetSpineHitList(), "spineHitListV", VIOLET);
    PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &protoShowerW.GetSpineHitList(), "spineHitListW", VIOLET);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &protoShowerU.GetShowerCore().GetStartPosition(), "showerStartU", BLUE, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &protoShowerV.GetShowerCore().GetStartPosition(), "showerStartV", BLUE, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &protoShowerW.GetShowerCore().GetStartPosition(), "showerStartW", BLUE, 2);
    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    //////////////////////////
    if (this->AreShowerStartsConsistent(protoShowerU, protoShowerV, protoShowerW))
    {
        consistency = Consistency::POSITION;
    }
    else if (this->AreDirectionsConsistent(protoShowerU, protoShowerV, protoShowerW))
    {
        consistency = Consistency::DIRECTION;
    }
    else
    {
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ProtoShowerMatchingTool::AreShowerStartsConsistent(const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW) const
{
    //std::cout << "AreShowerStartsConsistent? " << std::endl;

    const CartesianVector &showerStartU(protoShowerU.GetShowerCore().GetStartPosition());
    const CartesianVector &showerStartV(protoShowerV.GetShowerCore().GetStartPosition());
    const CartesianVector &showerStartW(protoShowerW.GetShowerCore().GetStartPosition());

    const float xSeparationUV(std::fabs(showerStartU.GetX() - showerStartV.GetX()));
    const float xSeparationUW(std::fabs(showerStartU.GetX() - showerStartW.GetX()));
    const float xSeparationVW(std::fabs(showerStartV.GetX() - showerStartW.GetX()));

    if ((xSeparationUV > m_maxXSeparation) || (xSeparationUW > m_maxXSeparation) || (xSeparationVW > m_maxXSeparation))
    {
        //std::cout << "X separation is bad.." << std::endl;
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    float chi2(0.f);
    CartesianVector projectionU(0.f, 0.f, 0.f), projectionV(0.f, 0.f, 0.f), projectionW(0.f, 0.f, 0.f);

    LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, showerStartV, showerStartW, projectionU, chi2);
    LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, showerStartW, showerStartU, projectionV, chi2);
    LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, showerStartU, showerStartV, projectionW, chi2);

    const float separationU((projectionU - showerStartU).GetMagnitude());
    const float separationV((projectionV - showerStartV).GetMagnitude());
    const float separationW((projectionW - showerStartW).GetMagnitude());

    ////////////////////////////
    /*
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &projectionU, "projectionShowerStartU", BLUE, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &projectionV, "projectionShowerStartV", BLUE, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &projectionW, "projectionShowerStartW", BLUE, 2);

    std::cout << "separationU: " << separationU << std::endl;
    std::cout << "separationV: " << separationV << std::endl;
    std::cout << "separationW: " << separationW << std::endl;

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    ////////////////////////////

    const float metric((separationU + separationV + separationW) / 3.f);

    return (metric < m_maxSeparation);
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool ProtoShowerMatchingTool::AreDirectionsConsistent(const ProtoShower &protoShowerU, const ProtoShower &protoShowerV, const ProtoShower &protoShowerW) const
{
    //std::cout << "Are Initial directions consistent?" << std::endl;

    const CartesianVector &directionU1(protoShowerU.GetConnectionPathway().GetStartDirection());
    const CartesianVector &directionV1(protoShowerV.GetConnectionPathway().GetStartDirection());
    const CartesianVector &directionW1(protoShowerW.GetConnectionPathway().GetStartDirection());

    
    if (this->AreDirectionsConsistent(protoShowerU.GetConnectionPathway().GetStartPosition(), protoShowerV.GetConnectionPathway().GetStartPosition(), protoShowerW.GetConnectionPathway().GetStartPosition(),
    directionU1, directionV1, directionW1))
        //if (this->AreDirectionsConsistent(directionU1, directionV1, directionW1))
    {
        return true;
    }
    else
    {
        //std::cout << "are shower spine directions consistent?" << std::endl;

        const bool isDownstream(
            protoShowerW.GetShowerCore().GetStartPosition().GetZ() > protoShowerW.GetConnectionPathway().GetStartPosition().GetZ());

        CartesianPointVector spinePositionsU, spinePositionsV, spinePositionsW;

        for (const CaloHit *const pCaloHit : protoShowerU.GetSpineHitList())
            spinePositionsU.push_back(pCaloHit->GetPositionVector());

        for (const CaloHit *const pCaloHit : protoShowerV.GetSpineHitList())
            spinePositionsV.push_back(pCaloHit->GetPositionVector());

        for (const CaloHit *const pCaloHit : protoShowerW.GetSpineHitList())
            spinePositionsW.push_back(pCaloHit->GetPositionVector());

        const TwoDSlidingFitResult spineFitU(&spinePositionsU, m_spineSlidingFitWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_U));
        const TwoDSlidingFitResult spineFitV(&spinePositionsV, m_spineSlidingFitWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_V));
        const TwoDSlidingFitResult spineFitW(&spinePositionsW, m_spineSlidingFitWindow, LArGeometryHelper::GetWirePitch(this->GetPandora(), TPC_VIEW_W));

        const CartesianVector directionU2(isDownstream ? spineFitU.GetGlobalMinLayerDirection() : spineFitU.GetGlobalMaxLayerDirection() * (-1.f));
        const CartesianVector directionV2(isDownstream ? spineFitV.GetGlobalMinLayerDirection() : spineFitV.GetGlobalMaxLayerDirection() * (-1.f));
        const CartesianVector directionW2(isDownstream ? spineFitW.GetGlobalMinLayerDirection() : spineFitW.GetGlobalMaxLayerDirection() * (-1.f));

        return this->AreDirectionsConsistent(protoShowerU.GetConnectionPathway().GetStartPosition(), protoShowerV.GetConnectionPathway().GetStartPosition(), protoShowerW.GetConnectionPathway().GetStartPosition(), directionU2, directionV2, directionW2);
        //return (this->AreDirectionsConsistent(directionU2, directionV2, directionW2));
    }

    return false;
}


//------------------------------------------------------------------------------------------------------------------------------------------

bool ProtoShowerMatchingTool::AreDirectionsConsistent(CartesianVector nuVertexU, CartesianVector nuVertexV,CartesianVector nuVertexW,
    CartesianVector directionU, CartesianVector directionV, CartesianVector directionW) const
{
    const CartesianVector wireAxis(0.f, 0.f, 1.f);

    float wireDeviationU(directionU.GetOpeningAngle(wireAxis));
    wireDeviationU = std::min(wireDeviationU, static_cast<float>(M_PI - wireDeviationU));

    float wireDeviationV(directionV.GetOpeningAngle(wireAxis));
    wireDeviationV = std::min(wireDeviationV, static_cast<float>(M_PI - wireDeviationV));

    float wireDeviationW(directionW.GetOpeningAngle(wireAxis));
    wireDeviationW = std::min(wireDeviationW, static_cast<float>(M_PI - wireDeviationW));

    float radians((2.f * M_PI) / 180.f);
    bool isIsochronous((wireDeviationU < radians) && (wireDeviationV < radians) && (wireDeviationW < radians));
    /*
    std::cout << "wireDeviationU: " << (wireDeviationU * 180.f / M_PI) << std::endl;
    std::cout << "wireDeviationV: " << (wireDeviationV * 180.f / M_PI) << std::endl;
    std::cout << "wireDeviationW: " << (wireDeviationW * 180.f / M_PI) << std::endl;
    */
    if (isIsochronous)
    {
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return true;
    }

    if (directionU.GetX() * directionV.GetX() < 0.f)
    {
        //std::cout << "UV issue" << std::endl;
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    if (directionU.GetX() * directionW.GetX() < 0.f)
    {
        //std::cout << "UW issue" << std::endl;
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    if (directionV.GetX() * directionW.GetX() < 0.f)
    {
        //std::cout << "VW issue" << std::endl;
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    const CartesianVector xAxis(1.f, 0.f, 0.f);

    const float cosOpeningAngleU(std::fabs(directionU.GetCosOpeningAngle(xAxis)));
    const float cosOpeningAngleV(std::fabs(directionV.GetCosOpeningAngle(xAxis)));
    const float cosOpeningAngleW(std::fabs(directionW.GetCosOpeningAngle(xAxis)));

    /*
    std::cout << "cosOpeningAngleU: " << cosOpeningAngleU << std::endl;
    std::cout << "cosOpeningAngleV: " << cosOpeningAngleV << std::endl;
    std::cout << "cosOpeningAngleW: " << cosOpeningAngleW << std::endl;
    */

    const float xSeparation(5.f);

    const CartesianVector uPoint(nuVertexU + (directionU * (xSeparation / cosOpeningAngleU)));
    const CartesianVector vPoint(nuVertexV + (directionV * (xSeparation / cosOpeningAngleV)));
    const CartesianVector wPoint(nuVertexW + (directionW * (xSeparation / cosOpeningAngleW)));

    /*
    std::cout << "uPoint: " << uPoint << std::endl;
    std::cout << "vPoint: " << vPoint << std::endl;
    std::cout << "wPoint: " << wPoint << std::endl;

    std::cout << "directionU: " << directionU << std::endl;
    std::cout << "directionV: " << directionV << std::endl;
    std::cout << "directionW: " << directionW << std::endl;
    */

    float chiSquared(0.f);

    CartesianVector projectionPointU(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, vPoint, wPoint, projectionPointU, chiSquared);

    CartesianVector projectionPointV(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, wPoint, uPoint, projectionPointV, chiSquared);

    CartesianVector projectionPointW(0.f, 0.f, 0.f);
    LArGeometryHelper::MergeTwoPositions(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, uPoint, vPoint, projectionPointW, chiSquared);

    // Project U
    const CartesianVector projectionU((projectionPointU - nuVertexU).GetUnitVector());

    // Project V
    const CartesianVector projectionV((projectionPointV - nuVertexV).GetUnitVector());

    // Project W
    const CartesianVector projectionW((projectionPointW - nuVertexW).GetUnitVector());

    ///////////////////////////////////////
    /*
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &uPoint, "uPoint", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &vPoint, "vPoint", BLACK, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &wPoint, "wPoint", BLACK, 2);

    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &projectionPointU, "projectionPointU", GREEN, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &projectionPointV, "projectionPointV", GREEN, 2);
    PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &projectionPointW, "projectionPointW", GREEN, 2);

    const CartesianVector endU(nuVertexU + (directionU * 30.f));
    const CartesianVector endV(nuVertexV + (directionV * 30.f));
    const CartesianVector endW(nuVertexW + (directionW * 30.f));

    const CartesianVector endU2(nuVertexU + (projectionU * 30.f));
    const CartesianVector endV2(nuVertexV + (projectionV * 30.f));
    const CartesianVector endW2(nuVertexW + (projectionW * 30.f));

    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &nuVertexU, &endU, "Direction U", BLACK, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &nuVertexV, &endV, "Direction V", BLACK, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &nuVertexW, &endW, "Direction W", BLACK, 2, 2);

    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &nuVertexU, &endU2, "Direction U", GREEN, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &nuVertexV, &endV2, "Direction V", GREEN, 2, 2);
    PandoraMonitoringApi::AddLineToVisualization(this->GetPandora(), &nuVertexW, &endW2, "Direction W", GREEN, 2, 2);
    */
    ///////////////////////////////////////

    float openingAngleU(directionU.GetOpeningAngle(projectionU) * 180.f / M_PI);
    float openingAngleV(directionV.GetOpeningAngle(projectionV) * 180.f / M_PI);
    float openingAngleW(directionW.GetOpeningAngle(projectionW) * 180.f / M_PI);

    /////////////////////
    /*
    std::cout << "openingAngleU: " << openingAngleU << std::endl;
    std::cout << "openingAngleV: " << openingAngleV << std::endl;
    std::cout << "openingAngleW: " << openingAngleW << std::endl;

    PandoraMonitoringApi::ViewEvent(this->GetPandora());
    */
    /////////////////////  

    const float metric((openingAngleU + openingAngleV + openingAngleW) / 3.f);

    if (metric > m_maxAngularDeviation)
        return false;

    return true;
}



//------------------------------------------------------------------------------------------------------------------------------------------

bool ProtoShowerMatchingTool::AreDirectionsConsistent(CartesianVector directionU, CartesianVector directionV, CartesianVector directionW) const
{
    const CartesianVector wireAxis(0.f, 0.f, 1.f);

    float wireDeviationU(directionU.GetOpeningAngle(wireAxis));
    wireDeviationU = std::min(wireDeviationU, static_cast<float>(M_PI - wireDeviationU));

    float wireDeviationV(directionV.GetOpeningAngle(wireAxis));
    wireDeviationV = std::min(wireDeviationV, static_cast<float>(M_PI - wireDeviationV));

    float wireDeviationW(directionW.GetOpeningAngle(wireAxis));
    wireDeviationW = std::min(wireDeviationW, static_cast<float>(M_PI - wireDeviationW));

    float radians((2.f * M_PI) / 180.f);
    bool isIsochronous((wireDeviationU < radians) && (wireDeviationV < radians) && (wireDeviationW < radians));

    if (isIsochronous)
        return true;

    if (isIsochronous)
    {
        // Enforce consistency
        directionU = CartesianVector(std::fabs(directionU.GetX()), 0.f, directionU.GetZ());
        directionV = CartesianVector(std::fabs(directionV.GetX()), 0.f, directionV.GetZ());
        directionW = CartesianVector(std::fabs(directionW.GetX()), 0.f, directionW.GetZ());
    }

    std::cout << "isIsochronous: " << isIsochronous << std::endl;

    if (directionU.GetX() * directionV.GetX() < 0.f)
    {
        //std::cout << "UV issue" << std::endl;
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    if (directionU.GetX() * directionW.GetX() < 0.f)
    {
        //std::cout << "UW issue" << std::endl;
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    if (directionV.GetX() * directionW.GetX() < 0.f)
    {
        //std::cout << "VW issue" << std::endl;
        //PandoraMonitoringApi::ViewEvent(this->GetPandora());
        return false;
    }

    const CartesianVector projectionU(LArGeometryHelper::MergeTwoDirections(this->GetPandora(), TPC_VIEW_V, TPC_VIEW_W, directionV, directionW));
    const CartesianVector projectionV(LArGeometryHelper::MergeTwoDirections(this->GetPandora(), TPC_VIEW_W, TPC_VIEW_U, directionW, directionU));
    const CartesianVector projectionW(LArGeometryHelper::MergeTwoDirections(this->GetPandora(), TPC_VIEW_U, TPC_VIEW_V, directionU, directionV));

    float openingAngleU(directionU.GetOpeningAngle(projectionU) * 180.f / M_PI);
    float openingAngleV(directionV.GetOpeningAngle(projectionV) * 180.f / M_PI);
    float openingAngleW(directionW.GetOpeningAngle(projectionW) * 180.f / M_PI);

    if (isIsochronous)
    {
        openingAngleU = std::min(openingAngleU, 180.f - openingAngleU);
        openingAngleV = std::min(openingAngleV, 180.f - openingAngleV);
        openingAngleW = std::min(openingAngleW, 180.f - openingAngleW);
    }

    std::cout << "openingAngleU: " << openingAngleU << std::endl;
    std::cout << "openingAngleV: " << openingAngleV << std::endl;
    std::cout << "openingAngleW: " << openingAngleW << std::endl;

    //PandoraMonitoringApi::ViewEvent(this->GetPandora());

    if ((openingAngleU > m_maxAngularDeviation) || (openingAngleV > m_maxAngularDeviation) || (openingAngleW > m_maxAngularDeviation))
        return false;

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode ProtoShowerMatchingTool::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SpineSlidingFitWindow", m_spineSlidingFitWindow));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxXSeparation", m_maxXSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxSeparation", m_maxSeparation));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxAngularDeviation", m_maxAngularDeviation));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
