/**
 *  @file   larpandoracontent/LArTwoDReco/LArClusterSplitting/WireSplittingAlgorithm.cc
 *
 *  @brief  Implementation of the wire splitting algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include "larpandoracontent/LArTwoDReco/LArClusterSplitting/WireSplittingAlgorithm.h"

using namespace pandora;

namespace lar_content
{

WireSplittingAlgorithm::WireSplittingAlgorithm() :
    m_minClusterLength(10.f)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode WireSplittingAlgorithm::DivideCaloHits(const Cluster *const pCluster, CaloHitList &firstHitList, CaloHitList &secondHitList) const
{
    if (LArClusterHelper::GetLengthSquared(pCluster) < m_minClusterLength * m_minClusterLength)
        return STATUS_CODE_NOT_FOUND;

    CaloHitList caloHitList;
    LArClusterHelper::GetAllHits(pCluster, caloHitList);

    // Define some things
    typedef std::pair<float, float> WireExtent;
    typedef std::pair<std::pair<WireExtent, WireExtent>, CaloHitList> HitGroup;
    std::vector<HitGroup> hitGroups;
    CaloHitList usedHits;
    
    // Create initial group
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit *>(pCaloHit));

        if (!pLArCaloHit)
            continue;

        const unsigned int minWireIntersect1(pLArCaloHit->GetMinIntersectWire1());
        const unsigned int maxWireIntersect1(pLArCaloHit->GetMaxIntersectWire1());
        const unsigned int minWireIntersect2(pLArCaloHit->GetMinIntersectWire2());
        const unsigned int maxWireIntersect2(pLArCaloHit->GetMaxIntersectWire2());        

        hitGroups.emplace_back(std::make_pair(std::make_pair(WireExtent(minWireIntersect1, maxWireIntersect1), WireExtent(minWireIntersect2, maxWireIntersect2)), CaloHitList{pCaloHit}));
        usedHits.emplace_back(pCaloHit);
        break;
    }

    // Now split into hit groups
    bool found(true);

    while(found)
    {
        found = false;
        
        for (const CaloHit *const pCaloHit : caloHitList)
        {
            if (std::find(usedHits.begin(), usedHits.end(), pCaloHit) != usedHits.end())
                continue;

            const LArCaloHit *const pLArCaloHit(dynamic_cast<const LArCaloHit *>(pCaloHit));

            if (!pLArCaloHit)
                continue;

            const unsigned int minWireIntersect1(pLArCaloHit->GetMinIntersectWire1());
            const unsigned int maxWireIntersect1(pLArCaloHit->GetMaxIntersectWire1());
            const unsigned int minWireIntersect2(pLArCaloHit->GetMinIntersectWire2());
            const unsigned int maxWireIntersect2(pLArCaloHit->GetMaxIntersectWire2());

            HitGroup &hitGroup(hitGroups.back());
            
            const unsigned int groupMinWireIntersect1(hitGroup.first.first.first);
            const unsigned int groupMaxWireIntersect1(hitGroup.first.first.second);
            const unsigned int groupMinWireIntersect2(hitGroup.first.second.first);
            const unsigned int groupMaxWireIntersect2(hitGroup.first.second.second);

            // check to see if the hit fits into the group, if not skip
            if ((minWireIntersect1 < groupMinWireIntersect1) && (maxWireIntersect1 < groupMinWireIntersect1))
                continue;

            if ((minWireIntersect1 > groupMaxWireIntersect1) && (maxWireIntersect1 > groupMaxWireIntersect1))
                continue;

            if ((minWireIntersect2 < groupMinWireIntersect2) && (maxWireIntersect2 < groupMinWireIntersect2))
                continue;

            if ((minWireIntersect2 > groupMaxWireIntersect2) && (maxWireIntersect2 > groupMaxWireIntersect2))
                continue;

            // if we get here then the hit fits into the group, so add it and update the group extent
            hitGroup.second.emplace_back(pCaloHit);
            hitGroup.first.first.first = std::min(groupMinWireIntersect1, minWireIntersect1);
            hitGroup.first.first.second = std::max(groupMaxWireIntersect1, maxWireIntersect1);
            hitGroup.first.second.first = std::min(groupMinWireIntersect2, minWireIntersect2);
            hitGroup.first.second.second = std::max(groupMaxWireIntersect2, maxWireIntersect2);
            usedHits.emplace_back(pCaloHit);
            found = true;
        }
    }

    // Visualise the hit group
    ClusterList visualiseClusters({pCluster});
    CaloHitList visualiseHitList;
    for (const CaloHit *const pCaloHit : caloHitList)
    {
        if (std::find(usedHits.begin(), usedHits.end(), pCaloHit) != usedHits.end())
            continue;
        
        CartesianVector position(pCaloHit->GetPositionVector());
        PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &position, "hit group", RED, 2);
        PandoraMonitoringApi::VisualizeClusters(this->GetPandora(), &visualiseClusters, "Cluster", BLUE);
        PandoraMonitoringApi::Pause(this->GetPandora());
    }

    if (caloHitList.size() != usedHits.size())
    {
        std::cout << "Event: " << this->GetPandora().GetEvent() << std::endl;
    }
    
    return STATUS_CODE_NOT_FOUND;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode WireSplittingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinClusterLength", m_minClusterLength));

    return ClusterSplittingAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_content
