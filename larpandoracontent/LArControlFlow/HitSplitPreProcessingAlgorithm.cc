/**
 *  @file   larpandoracontent/LArControlFlow/HitSplitPreProcessingAlgorithm.cc
 *
 *  @brief  Implementation of the list preparation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/HitSplitPreProcessingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

using namespace pandora;

namespace lar_content
{

HitSplitPreProcessingAlgorithm::HitSplitPreProcessingAlgorithm() :
    m_mipEquivalentCut(std::numeric_limits<float>::epsilon()),
    m_minCellLengthScale(std::numeric_limits<float>::epsilon()),
    m_maxCellLengthScale(std::numeric_limits<float>::max()),
    m_searchRegion1D(0.1f),
    m_maxEventHits(std::numeric_limits<unsigned int>::max()),
    m_onlyAvailableCaloHits(true),
    m_inputCaloHitListName("Input")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitSplitPreProcessingAlgorithm::Reset()
{
    m_processedHits.clear();
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitSplitPreProcessingAlgorithm::Run()
{
    if (!this->GetPandora().GetSettings()->SingleHitTypeClusteringMode())
    {
        std::cout << "HitSplitPreProcessingAlgorithm: expect Pandora to be configured in SingleHitTypeClusteringMode." << std::endl;
        return STATUS_CODE_FAILURE;
    }

    try
    {
        this->ProcessCaloHits();
    }
    catch (const StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_OUT_OF_RANGE == statusCodeException.GetStatusCode())
        {
            std::cout << "HitSplitPreProcessingAlgorithm: Excessive number of hits in event, skipping the reconstruction" << std::endl;
            this->PopulateVoidCaloHitLists();
        }
    }

    if (!m_currentCaloHitListReplacement.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_currentCaloHitListReplacement))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "HitSplitPreProcessingAlgorithm: could not replace current calo hit list with list named: " << m_currentCaloHitListReplacement << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitSplitPreProcessingAlgorithm::ProcessCaloHits()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListName, pCaloHitList));

    if (pCaloHitList->empty())
        return;

    if (pCaloHitList->size() > m_maxEventHits)
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    CaloHitList selectedCaloHitListU, selectedCaloHitListV, selectedCaloHitListW;

    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (m_processedHits.count(pCaloHit))
            continue;

        (void) m_processedHits.insert(pCaloHit);

        if (m_onlyAvailableCaloHits && !PandoraContentApi::IsAvailable(*this, pCaloHit))
            continue;

        if (pCaloHit->GetMipEquivalentEnergy() < m_mipEquivalentCut)
            continue;

        if (pCaloHit->GetInputEnergy() < std::numeric_limits<float>::epsilon())
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "HitSplitPreProcessingAlgorithm: found a hit with zero energy, will remove it" << std::endl;

            continue;
        }

        if ((pCaloHit->GetCellLengthScale() < m_minCellLengthScale) || (pCaloHit->GetCellLengthScale() > m_maxCellLengthScale))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            {
                std::cout << "HitSplitPreProcessingAlgorithm: found a hit with extent " << pCaloHit->GetCellLengthScale()
                          << ", require (" << m_minCellLengthScale << " - " << m_maxCellLengthScale << "), will remove it" << std::endl;
            }

            continue;
        }

        if (TPC_VIEW_U == pCaloHit->GetHitType())
        {
            selectedCaloHitListU.push_back(pCaloHit);
        }
        else if (TPC_VIEW_V == pCaloHit->GetHitType())
        {
            selectedCaloHitListV.push_back(pCaloHit);
        }
        else if (TPC_VIEW_W == pCaloHit->GetHitType())
        {
            selectedCaloHitListW.push_back(pCaloHit);
        }
    }

    CaloHitList constituentHitListU, constituentHitListV, constituentHitListW;
    SplitCaloHitList(selectedCaloHitListU, constituentHitListU);
    SplitCaloHitList(selectedCaloHitListV, constituentHitListV);
    SplitCaloHitList(selectedCaloHitListW, constituentHitListW);

    CaloHitList filteredCaloHitListU, filteredCaloHitListV, filteredCaloHitListW;
    this->GetFilteredCaloHitList(constituentHitListU, filteredCaloHitListU);
    this->GetFilteredCaloHitList(constituentHitListV, filteredCaloHitListV);
    this->GetFilteredCaloHitList(constituentHitListW, filteredCaloHitListW);

    CaloHitList filteredInputList;
    filteredInputList.insert(filteredInputList.end(), filteredCaloHitListU.begin(), filteredCaloHitListU.end());
    filteredInputList.insert(filteredInputList.end(), filteredCaloHitListV.begin(), filteredCaloHitListV.end());
    filteredInputList.insert(filteredInputList.end(), filteredCaloHitListW.begin(), filteredCaloHitListW.end());

    if (!filteredInputList.empty() && !m_filteredCaloHitListName.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredInputList, m_filteredCaloHitListName));

    if (!filteredCaloHitListU.empty() && !m_outputCaloHitListNameU.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredCaloHitListU, m_outputCaloHitListNameU));

    if (!filteredCaloHitListV.empty() && !m_outputCaloHitListNameV.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredCaloHitListV, m_outputCaloHitListNameV));

    if (!filteredCaloHitListW.empty() && !m_outputCaloHitListNameW.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredCaloHitListW, m_outputCaloHitListNameW));
}

//------------------------------------------------------------------------------------------------------------------------------------------


  void HitSplitPreProcessingAlgorithm::SplitCaloHitList(const CaloHitList &originalCaloHitList, CaloHitList &constituentCaloHitList)
{

    //input list will contain created calo hits
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListName, pCaloHitList));

    const CaloHitList initialInputCaloHitList(pCaloHitList->begin(), pCaloHitList->end());

    float hitWidth(0.3);

    for(const CaloHit *const pCaloHit : originalCaloHitList)
    {

      const unsigned int numberOfConstituents(floor(pCaloHit->GetCellSize1()/hitWidth) + 1);

      float constituentInputEnergy(pCaloHit->GetInputEnergy()/numberOfConstituents);
      float constituentMipEquivalentEnergy(pCaloHit->GetMipEquivalentEnergy()/numberOfConstituents);
      float constituentElectromagneticEnergy(pCaloHit->GetElectromagneticEnergy()/numberOfConstituents);
      float constituentHadronicEnergy(pCaloHit->GetHadronicEnergy()/numberOfConstituents);

      object_creation::CaloHit::Parameters caloHitParameters;
      caloHitParameters.m_positionVector = CartesianVector(0,0,0);
      caloHitParameters.m_expectedDirection = pCaloHit->GetExpectedDirection();
      caloHitParameters.m_cellNormalVector = pCaloHit->GetCellNormalVector();
      caloHitParameters.m_cellGeometry = pCaloHit->GetCellGeometry();
      caloHitParameters.m_cellSize0 = pCaloHit->GetCellSize0();
      caloHitParameters.m_cellSize1 = hitWidth;
      caloHitParameters.m_cellThickness = pCaloHit->GetCellThickness();
      caloHitParameters.m_nCellRadiationLengths = pCaloHit->GetNCellRadiationLengths();
      caloHitParameters.m_nCellInteractionLengths = pCaloHit->GetNCellInteractionLengths();
      caloHitParameters.m_time = pCaloHit->GetTime();
      caloHitParameters.m_inputEnergy = constituentInputEnergy;
      caloHitParameters.m_mipEquivalentEnergy = constituentMipEquivalentEnergy;
      caloHitParameters.m_electromagneticEnergy = constituentElectromagneticEnergy;
      caloHitParameters.m_hadronicEnergy = constituentHadronicEnergy;
      caloHitParameters.m_isDigital = pCaloHit->IsDigital();
      caloHitParameters.m_hitType = pCaloHit->GetHitType();
      caloHitParameters.m_hitRegion = pCaloHit->GetHitRegion();
      caloHitParameters.m_layer = pCaloHit->GetLayer();
      caloHitParameters.m_isInOuterSamplingLayer = pCaloHit->IsInOuterSamplingLayer();
      caloHitParameters.m_pParentAddress = pCaloHit->GetParentAddress();     

      CartesianVector positionAlongHit(pCaloHit->GetPositionVector());

      if(numberOfConstituents % 2 == 1)
      {
          caloHitParameters.m_positionVector = positionAlongHit;
          PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, object_creation::CaloHit::Create(this->GetPandora(), caloHitParameters));

	  for(unsigned int i(0); i < ((numberOfConstituents - 1)/2); ++i)
	  {
	      positionAlongHit += CartesianVector(hitWidth,0,0);
              caloHitParameters.m_positionVector = positionAlongHit;
              PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, object_creation::CaloHit::Create(this->GetPandora(), caloHitParameters));
	  }

	  positionAlongHit = pCaloHit->GetPositionVector();
	  for(unsigned int i(0); i < ((numberOfConstituents - 1)/2); ++i)
	  {
	      positionAlongHit -= CartesianVector(hitWidth,0,0);
              caloHitParameters.m_positionVector = positionAlongHit;
              PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, object_creation::CaloHit::Create(this->GetPandora(), caloHitParameters));
	  }
      }
      else
      {
          for(unsigned int i(0); i < (numberOfConstituents/2); ++i)
	  {
	      i == 0 ? positionAlongHit += CartesianVector(hitWidth/2,0,0) : positionAlongHit += CartesianVector(hitWidth,0,0);
              caloHitParameters.m_positionVector = positionAlongHit;
              PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, object_creation::CaloHit::Create(this->GetPandora(), caloHitParameters));
	  }
	  positionAlongHit = pCaloHit->GetPositionVector();
	  for(unsigned int i(0); i < (numberOfConstituents/2); ++i)
	  {
	      i == 0 ? positionAlongHit -= CartesianVector(hitWidth/2,0,0) :  positionAlongHit -= CartesianVector(hitWidth,0,0);
              caloHitParameters.m_positionVector = positionAlongHit;
              PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, object_creation::CaloHit::Create(this->GetPandora(), caloHitParameters));
	  }
      }

    }

    for(const CaloHit *const pCaloHit : *pCaloHitList)
    {
      //if calo hit is not in original list then must be a newly created hit
      if(std::find(initialInputCaloHitList.begin(), initialInputCaloHitList.end(), pCaloHit) == initialInputCaloHitList.end())
      {
	constituentCaloHitList.push_back(pCaloHit);
      }
    }

}


//------------------------------------------------------------------------------------------------------------------------------------------

void HitSplitPreProcessingAlgorithm::PopulateVoidCaloHitLists() noexcept
{
    CaloHitList emptyList;

    try
    {
        if (!m_filteredCaloHitListName.empty())
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, emptyList, m_filteredCaloHitListName));

        if (!m_outputCaloHitListNameU.empty())
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, emptyList, m_outputCaloHitListNameU));

        if (!m_outputCaloHitListNameV.empty())
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, emptyList, m_outputCaloHitListNameV));

        if (!m_outputCaloHitListNameW.empty())
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, emptyList, m_outputCaloHitListNameW));
    }
    catch (...)
    {
        std::cout << "HitSplitPreProcessingAlgorithm: Unable to create void calo hit lists" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void HitSplitPreProcessingAlgorithm::GetFilteredCaloHitList(const CaloHitList &inputList, CaloHitList &outputList)
{
    HitKDTree2D kdTree;
    HitKDNode2DList hitKDNode2DList;

    KDTreeBox hitsBoundingRegion2D = fill_and_bound_2d_kd_tree(inputList, hitKDNode2DList);
    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);

    // Remove hits that are in the same physical location!
    for (const CaloHit *const pCaloHit1 : inputList)
    {
        bool isUnique(true);
        KDTreeBox searchRegionHits(build_2d_kd_search_region(pCaloHit1, m_searchRegion1D, m_searchRegion1D));

        HitKDNode2DList found;
        kdTree.search(searchRegionHits, found);

        for (const auto &hit : found)
        {
            const CaloHit *const pCaloHit2(hit.data);

            if (pCaloHit1 == pCaloHit2)
                continue;

            const float displacementSquared((pCaloHit2->GetPositionVector() - pCaloHit1->GetPositionVector()).GetMagnitudeSquared());

            if (displacementSquared < std::numeric_limits<float>::epsilon())
            {
                const float deltaMip(pCaloHit2->GetMipEquivalentEnergy() > pCaloHit1->GetMipEquivalentEnergy());

                if ((deltaMip > std::numeric_limits<float>::epsilon()) ||
                    ((std::fabs(deltaMip) < std::numeric_limits<float>::epsilon()) && (outputList.end() != std::find(outputList.begin(), outputList.end(), pCaloHit2))))
                {
                    isUnique = false;
                    break;
                }
            }
        }

        if (isUnique)
        {
            outputList.push_back(pCaloHit1);
        }
        else
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "HitSplitPreProcessingAlgorithm: found two hits in same location, will remove lowest pulse height" << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode HitSplitPreProcessingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MipEquivalentCut", m_mipEquivalentCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MinCellLengthScale", m_minCellLengthScale));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxCellLengthScale", m_maxCellLengthScale));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "SearchRegion1D", m_searchRegion1D));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "MaxEventHits", m_maxEventHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "OnlyAvailableCaloHits", m_onlyAvailableCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle,
        "InputCaloHitListName", m_inputCaloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputCaloHitListNameU", m_outputCaloHitListNameU));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputCaloHitListNameV", m_outputCaloHitListNameV));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "OutputCaloHitListNameW", m_outputCaloHitListNameW));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "FilteredCaloHitListName", m_filteredCaloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle,
        "CurrentCaloHitListReplacement", m_currentCaloHitListReplacement));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
