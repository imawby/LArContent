/**
 *  @file   larpandoracontent/LArCheating/CheatingShowerStartFinderAlgorithm.cc
 *
 *  @brief  Implementation of the cheating shower start finder
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArCheating/CheatingShowerStartFinderAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArMCParticleHelper.h"

#include "PandoraMonitoringApi.h"

using namespace pandora;

namespace lar_content
{

CheatingShowerStartFinderAlgorithm::CheatingShowerStartFinderAlgorithm()
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingShowerStartFinderAlgorithm::Run()
{
    // Get MCParticle list
    const MCParticleList *pMCParticleList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_mcParticleListName, pMCParticleList));

    // Get CaloHit list
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_caloHitListName, pCaloHitList));

    std::cout << "pMCParticleList->size(): " << pMCParticleList->size() << std::endl;
    std::cout << "pCaloHitList->size(): " << pCaloHitList->size() << std::endl;

    const MCParticle *pMCPrimary(nullptr);

    for (const MCParticle *const pMCParticle : *pMCParticleList)
    {
        if (LArMCParticleHelper::IsPrimary(pMCParticle))
            pMCPrimary = pMCParticle;
    }

    if (!pMCPrimary)
    {
        std::cout << "No primary MC found!" << std::endl;
        return STATUS_CODE_SUCCESS;
    }

    // need map, to check we can see the particles..
    LArMCParticleHelper::MCRelationMap mcToSelfMap;
    LArMCParticleHelper::CaloHitToMCMap caloHitToMCMap;
    LArMCParticleHelper::MCContributionMap mcToCaloHitListMap; 
    LArMCParticleHelper::GetMCToSelfMap(pMCParticleList, mcToSelfMap);
    LArMCParticleHelper::GetMCParticleToCaloHitMatches(pCaloHitList, mcToSelfMap, caloHitToMCMap, mcToCaloHitListMap);

    HierarchyMap hierarchyMap;
    this->FillHierarchyMap(pMCPrimary, 1, hierarchyMap);

    for (auto &entry : hierarchyMap)
    {
        std::cout << "Generation: " << entry.first << std::endl;

        for (const MCParticle *const pMCParticle : entry.second)
        {
            if (mcToCaloHitListMap.find(pMCParticle) == mcToCaloHitListMap.end())
                continue;

            // have to split these up to view them -.-
            const CaloHitList &caloHitList(mcToCaloHitListMap.at(pMCParticle));
            CaloHitList caloHitListU, caloHitListV, caloHitListW;

            for (const CaloHit *const pCaloHit : caloHitList)
            {
                const HitType hitType(pCaloHit->GetHitType());
                CaloHitList &hitList(hitType == TPC_VIEW_U ? caloHitListU : hitType == TPC_VIEW_V ? caloHitListV : caloHitListW);
                hitList.push_back(pCaloHit);
            }

            const CartesianVector &vertex(pMCParticle->GetVertex());

            if (caloHitList.size() == 0)
                continue;

            std::cout << "caloHitList.size(): " << caloHitList.size() << std::endl;
            PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &caloHitListU, "MCHitsU", BLACK);
            PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &caloHitListV, "MCHitsV", BLACK);
            PandoraMonitoringApi::VisualizeCaloHits(this->GetPandora(), &caloHitListW, "MCHitsW", BLACK);
            PandoraMonitoringApi::AddMarkerToVisualization(this->GetPandora(), &vertex, "ElectronShowerVertex", BLACK, 2);
        }

        PandoraMonitoringApi::ViewEvent(this->GetPandora());
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------


void CheatingShowerStartFinderAlgorithm::FillHierarchyMap(const MCParticle *const pParentMCParticle, const int generation, 
    HierarchyMap &hierarchyMap)
{
    hierarchyMap[generation].push_back(pParentMCParticle);

    for (const MCParticle *const pMCParticle : pParentMCParticle->GetDaughterList())
        this->FillHierarchyMap(pMCParticle, generation + 1, hierarchyMap);
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode CheatingShowerStartFinderAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "MCParticleListName", m_mcParticleListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CaloHitListName", m_caloHitListName));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
