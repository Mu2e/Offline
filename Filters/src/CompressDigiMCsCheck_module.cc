// Check that the digi compression coming out of CompressDigiMCs worked
//
// Andy Edmonds, August 2018

#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

using namespace std;

namespace mu2e {

  class CompressDigiMCsCheck : public art::EDAnalyzer {
    art::InputTag _oldStrawDigiMCTag;
    art::InputTag _newStrawDigiMCTag;

    SimParticleTimeOffset _oldTOff;
    SimParticleTimeOffset _newTOff;

  public:
    explicit CompressDigiMCsCheck(const fhicl::ParameterSet& pset);
    virtual void analyze(const art::Event& event);
  };

  //================================================================
  CompressDigiMCsCheck::CompressDigiMCsCheck(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , _oldStrawDigiMCTag(pset.get<art::InputTag>("oldStrawDigiMCTag"))
    , _newStrawDigiMCTag(pset.get<art::InputTag>("newStrawDigiMCTag"))
    , _oldTOff(pset.get<fhicl::ParameterSet>("OldTimeOffsets", {}))
    , _newTOff(pset.get<fhicl::ParameterSet>("NewTimeOffsets", {}))
  {  }

  //================================================================
  void CompressDigiMCsCheck::analyze(const art::Event& event) {
    art::Handle<StrawDigiMCCollection> oldStrawDigiMCHandle;
    event.getByLabel(_oldStrawDigiMCTag, oldStrawDigiMCHandle);

    art::Handle<StrawDigiMCCollection> newStrawDigiMCHandle;
    event.getByLabel(_newStrawDigiMCTag, newStrawDigiMCHandle);

    if (oldStrawDigiMCHandle.isValid() && newStrawDigiMCHandle.isValid()) {
      const auto& oldStrawDigiMCs = *oldStrawDigiMCHandle;
      const auto& newStrawDigiMCs = *newStrawDigiMCHandle;

      _oldTOff.updateMap(event);
      _newTOff.updateMap(event);

      unsigned int n_old_straw_digi_mcs = oldStrawDigiMCs.size();
      unsigned int n_new_straw_digi_mcs = newStrawDigiMCs.size();

      if (n_old_straw_digi_mcs != n_new_straw_digi_mcs) {
	throw cet::exception("CompressDigiMCsCheck") << "Number of old and new StrawDigiMCs do not match" << std::endl;
      }
      
      for (unsigned int i_digi_mc = 0; i_digi_mc < n_old_straw_digi_mcs; ++i_digi_mc) {
	const auto& i_oldStrawDigiMC = oldStrawDigiMCs.at(i_digi_mc);
	const auto& i_newStrawDigiMC = newStrawDigiMCs.at(i_digi_mc);
	
	const auto& i_oldStepPointMC = *(i_oldStrawDigiMC.stepPointMC(StrawEnd::hv));
	const auto& i_newStepPointMC = *(i_newStrawDigiMC.stepPointMC(StrawEnd::hv));
	
	const auto& i_old_digi_mc_strawId = i_oldStrawDigiMC.strawId();
	const auto& i_new_digi_mc_strawId = i_newStrawDigiMC.strawId();
	if (i_old_digi_mc_strawId != i_new_digi_mc_strawId) {
	  throw cet::exception("CompressDigiMCsCheck") << "Old and new StrawDigiMC's StrawIds do not match" << std::endl;
	}
	
	const auto& i_old_step_strawId = i_oldStepPointMC.strawId();
	const auto& i_new_step_strawId = i_newStepPointMC.strawId();
	if (i_old_step_strawId != i_new_step_strawId) {
	  throw cet::exception("CompressDigiMCsCheck") << "Old and new StrawDigiMC's StepPointMC's StrawIds do not match" << std::endl;
	}
	
	if (i_new_step_strawId != i_new_digi_mc_strawId) {
	  throw cet::exception("CompressDigiMCsCheck") << "New StrawDigiMC's StrawId is inconsistent with its StepPointMC's StrawId" << std::endl;
	}
	
	double old_time = _oldTOff.timeWithOffsetsApplied(i_oldStepPointMC);
	double new_time = _newTOff.timeWithOffsetsApplied(i_newStepPointMC);
	if (std::fabs(old_time - new_time) > 1e-5) {
	  throw cet::exception("CompressDigiMCsCheck") << "Old and new StepPointMC times with offsets applied do not match" << std::endl;
	}
      }
    }
  } // analyze(event)

  //================================================================


} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CompressDigiMCsCheck);
