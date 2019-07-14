// Check that the digi compression coming out of CompressDigiMCs worked
//
// Andy Edmonds, August 2018

#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"

#include "MCDataProducts/inc/CaloShowerStepCollection.hh"

#include "TH1.h"

using namespace std;

namespace mu2e {

  class SeparateCaloShowerStepsCheck : public art::EDAnalyzer {
    art::InputTag _crystalStepsTag;
    art::InputTag _sipmStepsTag;

    TH1F* _histCrystalZ;
    TH1F* _histSiPMZ;

  public:
    explicit SeparateCaloShowerStepsCheck(const fhicl::ParameterSet& pset);
    virtual void analyze(const art::Event& event);
    virtual void beginJob();
  };

  //================================================================
  SeparateCaloShowerStepsCheck::SeparateCaloShowerStepsCheck(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , _crystalStepsTag(pset.get<art::InputTag>("crystalStepsTag"))
    , _sipmStepsTag(pset.get<art::InputTag>("sipmStepsTag"))
  {  }

  //===============================================================
  void SeparateCaloShowerStepsCheck::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;

    _histCrystalZ=tfs->make<TH1F>("histCrystalZ","Z-Position of \"calorimeter\" CaloShowerSteps", 5000,0,250);
    _histSiPMZ=tfs->make<TH1F>("histSiPMZ","Z-Position of \"calorimeterRO\" CaloShowerSteps", 5000,0,250);
  }

  //================================================================
  void SeparateCaloShowerStepsCheck::analyze(const art::Event& event) {
    art::Handle<CaloShowerStepCollection> crystalStepsHandle;
    event.getByLabel(_crystalStepsTag, crystalStepsHandle);

    if (crystalStepsHandle.isValid()) {
      for (const auto& i_crystalStep : *crystalStepsHandle) {
	_histCrystalZ->Fill(i_crystalStep.position().z());
      }
    }

    art::Handle<CaloShowerStepCollection> sipmStepsHandle;
    event.getByLabel(_sipmStepsTag, sipmStepsHandle);

    if (sipmStepsHandle.isValid()) {
      for (const auto& i_sipmStep : *sipmStepsHandle) {
	_histSiPMZ->Fill(i_sipmStep.position().z());
      }
    }

  } // analyze(event)

  //================================================================


} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SeparateCaloShowerStepsCheck);
