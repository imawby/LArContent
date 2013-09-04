/**
 *  @file   LArContent/src/LArClusterAssociation/LongitudinalExtensionAlgorithm.cc
 * 
 *  @brief  Implementation of the cluster association algorithm class.
 * 
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "LArClusterAssociation/LongitudinalExtensionAlgorithm.h"

#include "LArHelpers/LArVertexHelper.h"
#include "LArHelpers/LArClusterHelper.h"

using namespace pandora;

namespace lar
{

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::GetListOfCleanClusters( const ClusterList *const pClusterList, ClusterVector &clusterVector ) const
{
    for ( ClusterList::const_iterator iter = pClusterList->begin(), iterEnd = pClusterList->end(); iter != iterEnd; ++iter ) 
    {
        Cluster *pCluster = *iter;

        if ( 1 + pCluster->GetOuterPseudoLayer() - pCluster->GetInnerPseudoLayer() < 15 ) continue;

        const CartesianVector innerVertex = pCluster->GetCentroid( pCluster->GetInnerPseudoLayer() );
        const CartesianVector outerVertex = pCluster->GetCentroid( pCluster->GetOuterPseudoLayer() );

        if ( (outerVertex-innerVertex).GetMagnitudeSquared() < 25.0 ) continue;

        if ( LArClusterHelper::GetLayerOccupancy(pCluster) < 0.75 ) continue;

        clusterVector.push_back(pCluster);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillAssociationMatrix( const ClusterVector& clusterVector, ClusterAssociationMatrix& clusterAssociationMatrix) const
{
    // Generate a list of clean pointing clusters
    LArPointingClusterMap pointingClusterMap;
    LArPointingClusterList pointingClusterList;

    for (ClusterVector::const_iterator iter = clusterVector.begin(), iterEnd = clusterVector.end(); iter != iterEnd; ++iter)
    {
        pointingClusterMap.insert( std::pair<Cluster*,LArPointingCluster>(*iter,LArPointingCluster(*iter)) );
    }

    for (LArPointingClusterMap::const_iterator iter = pointingClusterMap.begin(), iterEnd = pointingClusterMap.end(); iter != iterEnd; ++iter)
    {
        pointingClusterList.push_back(iter->second);
    }


    // Form the matrix of associations
    for (LArPointingClusterList::const_iterator iterI = pointingClusterList.begin(), iterEndI = pointingClusterList.end(); iterI != iterEndI; ++iterI)
    {
        const LArPointingCluster& clusterI = *iterI;

        for (LArPointingClusterList::const_iterator iterJ = iterI, iterEndJ = pointingClusterList.end(); iterJ != iterEndJ; ++iterJ)
        {
            const LArPointingCluster& clusterJ = *iterJ;

            if (clusterI.GetCluster() == clusterJ.GetCluster())
                continue;

            this->FillAssociationMatrix(clusterI, clusterJ, clusterAssociationMatrix);
	}
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillAssociationMatrix( const LArPointingCluster& clusterI, const LArPointingCluster& clusterJ, ClusterAssociationMatrix& clusterAssociationMatrix ) const
{
    // Sanity Check
    if ( clusterI.GetCluster() == clusterJ.GetCluster() ) return ;

    // Check association of closest two end points
    const LArPointingCluster::Vertex& innerVertexI = clusterI.GetInnerVertex();
    const LArPointingCluster::Vertex& outerVertexI = clusterI.GetOuterVertex();
    const LArPointingCluster::Vertex& innerVertexJ = clusterJ.GetInnerVertex();
    const LArPointingCluster::Vertex& outerVertexJ = clusterJ.GetOuterVertex();

    const float distSquared_innerI_to_innerJ = (innerVertexI.GetPosition()-innerVertexJ.GetPosition()).GetMagnitudeSquared();
    const float distSquared_innerI_to_outerJ = (innerVertexI.GetPosition()-outerVertexJ.GetPosition()).GetMagnitudeSquared();
    const float distSquared_outerI_to_innerJ = (outerVertexI.GetPosition()-innerVertexJ.GetPosition()).GetMagnitudeSquared();
    const float distSquared_outerI_to_outerJ = (outerVertexI.GetPosition()-outerVertexJ.GetPosition()).GetMagnitudeSquared();

    // (a) Check association of inner I <-> inner J
    if ( distSquared_innerI_to_innerJ < std::min( distSquared_outerI_to_outerJ ,std::min( distSquared_innerI_to_outerJ, distSquared_outerI_to_innerJ ) )   
      && distSquared_outerI_to_outerJ > std::max( distSquared_innerI_to_innerJ, std::max( distSquared_innerI_to_outerJ, distSquared_outerI_to_innerJ ) ) )
    {
        return this->FillAssociationMatrix( innerVertexI, innerVertexJ, clusterAssociationMatrix  );
    }
    
    // (b) Check association of inner I <-> outer J
    if ( distSquared_innerI_to_outerJ < std::min( distSquared_outerI_to_innerJ, std::min( distSquared_outerI_to_outerJ, distSquared_innerI_to_innerJ ) )  
      && distSquared_outerI_to_innerJ > std::max( distSquared_innerI_to_outerJ, std::max( distSquared_outerI_to_outerJ, distSquared_innerI_to_innerJ ) ) )
    {
        return this->FillAssociationMatrix( innerVertexI, outerVertexJ, clusterAssociationMatrix );
    }
    
    // (c) Check association of outer I <-> inner J
    if ( distSquared_outerI_to_innerJ < std::min( distSquared_innerI_to_outerJ, std::min( distSquared_innerI_to_innerJ, distSquared_outerI_to_outerJ ) ) 
      && distSquared_innerI_to_outerJ > std::max( distSquared_outerI_to_innerJ, std::max( distSquared_innerI_to_innerJ, distSquared_outerI_to_outerJ ) ) ) 
    {
        return this->FillAssociationMatrix( outerVertexI, innerVertexJ, clusterAssociationMatrix );
    }
    
    // (d) Check association of outer I <-> outer J
    if ( distSquared_outerI_to_outerJ < std::min( distSquared_innerI_to_innerJ, std::min( distSquared_innerI_to_outerJ, distSquared_outerI_to_innerJ ) ) 
      && distSquared_innerI_to_innerJ > std::max( distSquared_outerI_to_outerJ, std::max( distSquared_innerI_to_outerJ, distSquared_outerI_to_innerJ ) ) )
    {
        return this->FillAssociationMatrix( outerVertexI, outerVertexJ, clusterAssociationMatrix );
    }

    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LongitudinalExtensionAlgorithm::FillAssociationMatrix( const LArPointingCluster::Vertex &clusterVertexI, const LArPointingCluster::Vertex &clusterVertexJ, ClusterAssociationMatrix& clusterAssociationMatrix ) const
{  
    // Sanity Check
    if ( clusterVertexI.GetCluster() == clusterVertexJ.GetCluster() ) return;

   
    // Check that linear fits have a reasonable RMS
    const float clusterRmsI = clusterVertexI.GetRms();
    const float clusterRmsJ = clusterVertexJ.GetRms();

    if( clusterRmsI > 0.5 || clusterRmsJ > 0.5 ) return;


    // Check that new layer occupancy would be reasonable
    const Cluster* pClusterI = clusterVertexI.GetCluster();
    const Cluster* pClusterJ = clusterVertexJ.GetCluster();

    const float clusterOccupancyI = LArClusterHelper::GetLayerOccupancy(pClusterI);
    const float clusterOccupancyJ = LArClusterHelper::GetLayerOccupancy(pClusterJ);

    if ( LArClusterHelper::GetLayerOccupancy(pClusterI, pClusterJ) < 0.75 ) return;


    // Calculate lengths of clusters
    const CartesianVector& innerCentroidI = pClusterI->GetCentroid( pClusterI->GetInnerPseudoLayer() );
    const CartesianVector& outerCentroidI = pClusterI->GetCentroid( pClusterI->GetOuterPseudoLayer() );

    const CartesianVector& innerCentroidJ = pClusterJ->GetCentroid( pClusterJ->GetInnerPseudoLayer() );
    const CartesianVector& outerCentroidJ = pClusterJ->GetCentroid( pClusterJ->GetOuterPseudoLayer() );

    const float clusterLengthI = clusterOccupancyI * (outerCentroidI - innerCentroidI).GetMagnitude();
    const float clusterLengthJ = clusterOccupancyJ * (outerCentroidJ - innerCentroidJ).GetMagnitude();
  

    // Association types
    ClusterAssociation::AssociationType associationType(ClusterAssociation::NOTHING);
    ClusterAssociation::StrengthType strengthType(ClusterAssociation::UNASSOCIATED);

    
    // Spatial resolution, hard-coded for now
    const float m_spatialResolution = 1.16; // cm


    // Requirements on cluster end points
    const CartesianVector& vertexPositionI = clusterVertexI.GetPosition();
    const CartesianVector& vertexPositionJ = clusterVertexJ.GetPosition();

    const CartesianVector& vertexDirectionI = clusterVertexI.GetDirection();
    const CartesianVector& vertexDirectionJ = clusterVertexJ.GetDirection();

    float distanceSquared = (vertexPositionI-vertexPositionJ).GetMagnitudeSquared();
    float cosTheta = -vertexDirectionI.GetDotProduct(vertexDirectionJ);
    
    if ( distanceSquared < 2.0 * m_spatialResolution * m_spatialResolution ) 
    {  
        associationType = ClusterAssociation::NODE;
        strengthType = ClusterAssociation::WEAK;
    
        if ( distanceSquared < m_spatialResolution * m_spatialResolution ) 
        { 
	    if ( cosTheta > 0.906 )
	    {
                strengthType = ClusterAssociation::STANDARD;

                if ( cosTheta > 0.966 ) strengthType = ClusterAssociation::STRONG;
	    }
	}
    }
   

    // Requirements on impact parameters  
    if ( clusterLengthI > 10.f && clusterLengthJ > 10.f )
    {
        float rT1(0.f), rL1(0.f), rT2(0.f), rL2(0.f);  

        // calculate impact parameters
        if ( strengthType < ClusterAssociation::STRONG )
        {
            const float cosThetaI = (vertexPositionI-vertexPositionJ).GetUnitVector().GetDotProduct(vertexDirectionI);
            const float cosThetaJ = (vertexPositionJ-vertexPositionI).GetUnitVector().GetDotProduct(vertexDirectionJ);

            LArVertexHelper::GetImpactParameters(vertexPositionI, vertexDirectionI, vertexPositionJ, rL1, rT1);
            LArVertexHelper::GetImpactParameters(vertexPositionJ, vertexDirectionJ, vertexPositionI, rL2, rT2);

            if ( cosTheta > 0.985 && std::fabs(cosThetaI) > 0.25 && std::fabs(cosThetaJ) > 0.25 
              && std::fabs(rL1) < 10.f && std::fabs(rL2) < 10.f && rT1 < 2.0 * m_spatialResolution && rT2 < 2.0 * m_spatialResolution )
            {
                associationType = ClusterAssociation::EMISSION;
                strengthType = ClusterAssociation::STRONG;
	    }
	}


        // remove 10 layers and re-calculate impact parameters
        if ( strengthType < ClusterAssociation::STRONG )
        {
            float truncatedRmsI(0.f);
            CartesianVector truncatedPositionI(0.f,0.f,0.f);
            CartesianVector truncatedDirectionI(0.f,0.f,0.f);    
    
            clusterVertexI.GetVertex(10, truncatedPositionI, truncatedDirectionI, truncatedRmsI);

            cosTheta = -truncatedDirectionI.GetDotProduct(vertexDirectionJ);

            LArVertexHelper::GetImpactParameters(truncatedPositionI, truncatedDirectionI, vertexPositionJ, rL1, rT1);
            LArVertexHelper::GetImpactParameters(vertexPositionJ, vertexDirectionJ, truncatedPositionI, rL2, rT2);

            if ( cosTheta > 0.985 && truncatedRmsI < 0.5
              && std::fabs(rL1) < 10.f && std::fabs(rL2) < 10.f && rT1 < m_spatialResolution && rT2 < m_spatialResolution )
            {
                associationType = ClusterAssociation::EMISSION;
                strengthType = ClusterAssociation::STRONG;
	    }
        }

        if ( strengthType < ClusterAssociation::STRONG )
	{
            float truncatedRmsJ(0.f);
            CartesianVector truncatedPositionJ(0.f,0.f,0.f);
            CartesianVector truncatedDirectionJ(0.f,0.f,0.f);    

            clusterVertexJ.GetVertex(10, truncatedPositionJ, truncatedDirectionJ, truncatedRmsJ);

            cosTheta = -truncatedDirectionJ.GetDotProduct(vertexDirectionI);

            LArVertexHelper::GetImpactParameters(truncatedPositionJ, truncatedDirectionJ, vertexPositionI, rL1, rT1);
            LArVertexHelper::GetImpactParameters(vertexPositionI, vertexDirectionI, truncatedPositionJ, rL2, rT2);

            if ( cosTheta > 0.985 && truncatedRmsJ < 0.5
             && std::fabs(rL1) < 10.f && std::fabs(rL2) < 10.f && rT1 < m_spatialResolution && rT2 < m_spatialResolution )
            { 
                associationType = ClusterAssociation::EMISSION;
                strengthType = ClusterAssociation::STRONG;
	    }
	}
    }

    // make entry in association matrix
    if ( strengthType > ClusterAssociation::UNASSOCIATED ) 
    {
        ClusterAssociation::VertexType vertexTypeI(ClusterAssociation::NONE);
        ClusterAssociation::VertexType vertexTypeJ(ClusterAssociation::NONE);

        if ( clusterVertexI.IsInner()==true ) vertexTypeI = ClusterAssociation::INNER; else vertexTypeI = ClusterAssociation::OUTER;
        if ( clusterVertexJ.IsInner()==true ) vertexTypeJ = ClusterAssociation::INNER; else vertexTypeJ = ClusterAssociation::OUTER;

        Cluster* pClusterI = clusterVertexI.GetCluster();
        Cluster* pClusterJ = clusterVertexJ.GetCluster();

        clusterAssociationMatrix[pClusterI][pClusterJ] = ClusterAssociation( vertexTypeI, vertexTypeJ, associationType, strengthType, clusterLengthJ );
        clusterAssociationMatrix[pClusterJ][pClusterI] = ClusterAssociation( vertexTypeJ, vertexTypeI, associationType, strengthType, clusterLengthI );

// ClusterList tempListI, tempListJ;
// tempListI.insert(pClusterI);
// tempListJ.insert(pClusterJ);
// CartesianVector positionI = clusterVertexI.GetPosition();
// CartesianVector positionJ = clusterVertexJ.GetPosition();
// PandoraMonitoringApi::SetEveDisplayParameters(0, 0, -1.f, 1.f);
// PandoraMonitoringApi::VisualizeClusters(&tempListI, "ClusterI",  GREEN);
// PandoraMonitoringApi::VisualizeClusters(&tempListJ, "ClusterJ",  BLUE);
// PandoraMonitoringApi::AddMarkerToVisualization(&positionI,"VertexI",RED, 1.);
// PandoraMonitoringApi::AddMarkerToVisualization(&positionJ,"VertexJ",RED, 1.);
// PandoraMonitoringApi::ViewEvent();
        
    }

}

//------------------------------------------------------------------------------------------------------------------------------------------

bool LongitudinalExtensionAlgorithm::AreClustersAssociated( Cluster* pCluster1, Cluster* pCluster2, ClusterAssociationMatrix& clusterAssociationMatrix ) const
{

    // Search for the pCluster1 <-> pCluster2 association
    ClusterAssociationMatrix::iterator iter1A = clusterAssociationMatrix.find(pCluster1);

    if ( clusterAssociationMatrix.end() == iter1A ) return false;

    ClusterAssociationMap& clusterAssociationMap1(iter1A->second);

    ClusterAssociationMap::iterator iter1B = clusterAssociationMap1.find(pCluster2);

    if ( clusterAssociationMap1.end() == iter1B ) return false;

    ClusterAssociation& clusterAssociationCandidate1(iter1B->second);

    if ( clusterAssociationCandidate1.GetStrength() < ClusterAssociation::STANDARD ) return false;


    // Bail out if there is a stronger association, or if this is a prong
    for ( ClusterAssociationMap::iterator iter1C = clusterAssociationMap1.begin(), iter1C_end = clusterAssociationMap1.end(); iter1C != iter1C_end; ++iter1C )
    {
        const Cluster* pClusterAlternative1(iter1C->first);

        if ( pClusterAlternative1 == pCluster2 ) continue;

        ClusterAssociation& clusterAssociationAlternative1(iter1C->second);

        if ( clusterAssociationCandidate1.GetParent() == clusterAssociationAlternative1.GetParent() )
        {
	    if ( clusterAssociationAlternative1.GetAssociation() == ClusterAssociation::NODE ) 
                return false;

	    if ( ( clusterAssociationAlternative1.GetStrength() > clusterAssociationCandidate1.GetStrength() )
              || ( clusterAssociationAlternative1.GetStrength() == clusterAssociationCandidate1.GetStrength()
		&& clusterAssociationAlternative1.GetFigureOfMerit() > clusterAssociationCandidate1.GetFigureOfMerit() ) )
	        return false;
        }
    }


    // Search for the pCluster2 <-> pCluster1 association
    ClusterAssociationMatrix::iterator iter2A = clusterAssociationMatrix.find(pCluster2);

    if ( clusterAssociationMatrix.end() == iter2A ) assert(0);

    ClusterAssociationMap& clusterAssociationMap2(iter2A->second);

    ClusterAssociationMap::iterator iter2B = clusterAssociationMap2.find(pCluster1);

    if ( clusterAssociationMap2.end() == iter2B ) assert(0);

    ClusterAssociation& clusterAssociationCandidate2(iter2B->second);

    if ( clusterAssociationCandidate2.GetStrength() < ClusterAssociation::STANDARD ) assert(0);


    // Bail out if there is a stronger association, or if this is a prong
    for ( ClusterAssociationMap::iterator iter2C = clusterAssociationMap2.begin(), iter2C_end = clusterAssociationMap2.end(); iter2C != iter2C_end; ++iter2C )
    {
        const Cluster* pClusterAlternative2(iter2C->first);

        if ( pClusterAlternative2 == pCluster1 ) continue;

        ClusterAssociation& clusterAssociationAlternative2(iter2C->second);

        if ( clusterAssociationCandidate2.GetParent() == clusterAssociationAlternative2.GetParent() )
        {
	    if ( clusterAssociationAlternative2.GetAssociation() == ClusterAssociation::NODE ) 
	        return false;

	    if ( ( clusterAssociationAlternative2.GetStrength() > clusterAssociationCandidate2.GetStrength() )
              || ( clusterAssociationAlternative2.GetStrength() == clusterAssociationCandidate2.GetStrength()
		&& clusterAssociationAlternative2.GetFigureOfMerit() > clusterAssociationCandidate2.GetFigureOfMerit() ) )
	        return false;
        }
    }
    
    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LongitudinalExtensionAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    
    return ClusterMergingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar
