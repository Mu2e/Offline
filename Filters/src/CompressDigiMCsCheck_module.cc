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

#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "DataProducts/inc/IndexMap.hh"
#include "MCDataProducts/inc/CrvDigiMC.hh"
#include "MCDataProducts/inc/CaloShowerSim.hh"

using namespace std;

namespace mu2e {

  class CompressDigiMCsCheck : public art::EDAnalyzer {
    art::InputTag _oldStrawDigiMCTag;
    art::InputTag _newStrawDigiMCTag;

    SimParticleTimeOffset _oldTOff;
    SimParticleTimeOffset _newTOff;

    art::InputTag _strawDigiMCIndexMapTag;
    IndexMap _strawDigiMCIndexMap;

    art::InputTag _oldCrvDigiMCTag;
    art::InputTag _newCrvDigiMCTag;
    art::InputTag _crvDigiMCIndexMapTag;
    IndexMap _crvDigiMCIndexMap;

    art::InputTag _oldCaloShowerSimTag;
    art::InputTag _newCaloShowerSimTag;

    bool _checkTrackerDuplicateSteps;
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
    , _strawDigiMCIndexMapTag(pset.get<art::InputTag>("strawDigiMCIndexMapTag", ""))
    , _oldCrvDigiMCTag(pset.get<art::InputTag>("oldCrvDigiMCTag"))
    , _newCrvDigiMCTag(pset.get<art::InputTag>("newCrvDigiMCTag"))
    , _crvDigiMCIndexMapTag(pset.get<art::InputTag>("crvDigiMCIndexMapTag", ""))
    , _oldCaloShowerSimTag(pset.get<art::InputTag>("oldCaloShowerSimTag"))
    , _newCaloShowerSimTag(pset.get<art::InputTag>("newCaloShowerSimTag"))
    , _checkTrackerDuplicateSteps(pset.get<bool>("checkTrackerDuplicateSteps",true))
  {  }

  //================================================================
  void CompressDigiMCsCheck::analyze(const art::Event& event) {

    ////////////////////////////////////
    // Check StrawDigiMCs
    art::Handle<StrawDigiMCCollection> oldStrawDigiMCHandle;
    event.getByLabel(_oldStrawDigiMCTag, oldStrawDigiMCHandle);

    art::Handle<StrawDigiMCCollection> newStrawDigiMCHandle;
    event.getByLabel(_newStrawDigiMCTag, newStrawDigiMCHandle);

    if (_strawDigiMCIndexMapTag != "") {
      art::Handle<IndexMap> strawDigiMCIndexMapHandle;
      event.getByLabel(_strawDigiMCIndexMapTag, strawDigiMCIndexMapHandle);
      _strawDigiMCIndexMap = *strawDigiMCIndexMapHandle;
    }

    if (oldStrawDigiMCHandle.isValid() && newStrawDigiMCHandle.isValid()) {
      const auto& oldStrawDigiMCs = *oldStrawDigiMCHandle;
      const auto& newStrawDigiMCs = *newStrawDigiMCHandle;

      _oldTOff.updateMap(event);
      _newTOff.updateMap(event);

      unsigned int n_old_straw_digi_mcs = oldStrawDigiMCs.size();
      unsigned int n_new_straw_digi_mcs = newStrawDigiMCs.size();

      // Only check for this if we haven't reduced the number of StrawDigiMCs
      if (_strawDigiMCIndexMapTag == "" && n_old_straw_digi_mcs != n_new_straw_digi_mcs) {
	throw cet::exception("CompressDigiMCsCheck") << "Number of old and new StrawDigiMCs do not match" << std::endl;
      }
      
      for (unsigned int i_old_digi_mc = 0; i_old_digi_mc < n_old_straw_digi_mcs; ++i_old_digi_mc) {
	const auto& i_oldStrawDigiMC = oldStrawDigiMCs.at(i_old_digi_mc);
	unsigned int i_new_digi_mc = i_old_digi_mc;
	if (_strawDigiMCIndexMapTag != "") {
	  if (_strawDigiMCIndexMap.checkInMap(i_old_digi_mc)) {
	    i_new_digi_mc = _strawDigiMCIndexMap.getCondensedIndex(i_old_digi_mc);
	  }
	  else {
	    continue; // to next old digi MC, since this one was compressed out...
	  }
	}
	const auto& i_newStrawDigiMC = newStrawDigiMCs.at(i_new_digi_mc);
	
	const auto& i_oldStepPointMC = i_oldStrawDigiMC.stepPointMC(StrawEnd::hv);
	if (!i_oldStepPointMC.isAvailable()) {
	  continue; // this is a null step point
	}
	const auto& i_newStepPointMC = i_newStrawDigiMC.stepPointMC(StrawEnd::hv);
	
	const auto& i_old_digi_mc_strawId = i_oldStrawDigiMC.strawId();
	const auto& i_new_digi_mc_strawId = i_newStrawDigiMC.strawId();
	if (i_old_digi_mc_strawId != i_new_digi_mc_strawId) {
	  throw cet::exception("CompressDigiMCsCheck") << "Old and new StrawDigiMC's StrawIds do not match" << std::endl;
	}
	
	const auto& i_old_step_strawId = i_oldStepPointMC->strawId();
	const auto& i_new_step_strawId = i_newStepPointMC->strawId();
	if (i_old_step_strawId != i_new_step_strawId) {
	  throw cet::exception("CompressDigiMCsCheck") << "Old and new StrawDigiMC's StepPointMC's StrawIds do not match" << std::endl;
	}
	
	if (i_new_step_strawId != i_new_digi_mc_strawId) {
	  throw cet::exception("CompressDigiMCsCheck") << "New StrawDigiMC's StrawId is inconsistent with its StepPointMC's StrawId" << std::endl;
	}
	
	double old_time = _oldTOff.timeWithOffsetsApplied(*i_oldStepPointMC);
	double new_time = _newTOff.timeWithOffsetsApplied(*i_newStepPointMC);
	if (std::fabs(old_time - new_time) > 1e-5) {
	  throw cet::exception("CompressDigiMCsCheck") << "Old and new StepPointMC times with offsets applied do not match (StrawDigiMC)" << std::endl;
	}

	// Check for tracker duplicate steps, there are two different cases. What should be happening is something like the following:
	// We have five StepPointMCs: A, B, C, D, E
	// StrawDigiMC has seven StepPtrs: HVStepPtr-->A   CalStepPtr-->A  WaveformStepPtrs-->A, B, C, D, E
	if (_checkTrackerDuplicateSteps) {
	  // Case 1. In this case we have accidentally treated the HVStepPtr and CalStepPtr as separate and duplicated step A (A' and A'')
	  //         so we end up with seven StepPointsMCs: A, B, C, D, E, A', A''
	  //         and the StepPtrs all point to different: HVStepPtr-->A'   CalStepPtr-->A''   WaveformStepPtrs-->A, B, C, D, E
	  // This will not trigger the exception in Case 2
	  // If we only access steps through digis then you will get the correct information
	  // You cannot loop over the StepPointMCCollection however
	  // Here we check that the HV and Cal StepPtrs also exist in the WaveformStepPtrs
	  const auto& i_newStepPointMCCal = i_newStrawDigiMC.stepPointMC(StrawEnd::cal);
	  bool identical_cal_ptr = false;
	  bool identical_hv_ptr = false;
	  for (const auto& i_triggerStepPointPtr : i_newStrawDigiMC.stepPointMCs()) {
	    if (i_triggerStepPointPtr == i_newStepPointMC) {
	      identical_hv_ptr = true;
	    }
	    if (i_triggerStepPointPtr == i_newStepPointMCCal) {
	      identical_cal_ptr = true;
	    }
	  }
	  if (! (identical_hv_ptr && identical_cal_ptr) ) {
	    throw cet::exception("CompressDigiMCsCheck") << "Trigger StepPointMCs in StrawDigiMCs are not identical to any StepPointMC in the waveform. This could indicate a duplication of StepPointMCs" << std::endl;
	  }

	  // Case 2: In this case, we have looped over an StepPointMCCollection with duplicated steps
	  //         (and have not duplicated things further)
	  //         StepPointMCs: A, B, C, D, E, A', A''
	  //         and the waveform will now include the duplicated steps:
	  //         StrawDigiMC: HVStepPtr-->A', CalStepPtr-->A'', WaveformStepPtrs-->A, B, C, D, E, A', A''
	  // This will not trigger the exception in Case 1
	  // Here we check that there are no identical waveform steps
	  for (const auto& i_stepPointPtr : i_newStrawDigiMC.stepPointMCs()) {
	    for (const auto& j_stepPointPtr : i_newStrawDigiMC.stepPointMCs()) {
	      if (i_stepPointPtr == j_stepPointPtr) {
		continue; // these will obviously match
	      }

	      if ( (i_stepPointPtr->volumeId() == j_stepPointPtr->volumeId()) &&
		   (i_stepPointPtr->totalEDep() == j_stepPointPtr->totalEDep())
		   ) {
		throw cet::exception("CompressDigiMCsCheck") << "Two StepPointMCs in StrawDigiMC waveform are identical" << std::endl;
	      }
	    }
	  }
	}
      }
    }

    //////////////////////////////////////
    // Check CrvDigiMCs
    art::Handle<CrvDigiMCCollection> oldCrvDigiMCHandle;
    event.getByLabel(_oldCrvDigiMCTag, oldCrvDigiMCHandle);

    art::Handle<CrvDigiMCCollection> newCrvDigiMCHandle;
    event.getByLabel(_newCrvDigiMCTag, newCrvDigiMCHandle);

    if (_crvDigiMCIndexMapTag != "") {
      art::Handle<IndexMap> crvDigiMCIndexMapHandle;
      event.getByLabel(_crvDigiMCIndexMapTag, crvDigiMCIndexMapHandle);
      _crvDigiMCIndexMap = *crvDigiMCIndexMapHandle;
    }

    if (oldCrvDigiMCHandle.isValid() && newCrvDigiMCHandle.isValid()) {
      const auto& oldCrvDigiMCs = *oldCrvDigiMCHandle;
      const auto& newCrvDigiMCs = *newCrvDigiMCHandle;

      _oldTOff.updateMap(event);
      _newTOff.updateMap(event);

      unsigned int n_old_crv_digi_mcs = oldCrvDigiMCs.size();
      unsigned int n_new_crv_digi_mcs = newCrvDigiMCs.size();

      // Only check for this if we haven't reduced the number of CrvDigiMCs (e.g. in digi production)
      if (_crvDigiMCIndexMapTag == "" && n_old_crv_digi_mcs != n_new_crv_digi_mcs) {
	throw cet::exception("CompressDigiMCsCheck") << "Number of old and new CrvDigiMCs do not match" << std::endl;
      }
      
      for (unsigned int i_old_digi_mc = 0; i_old_digi_mc < n_old_crv_digi_mcs; ++i_old_digi_mc) {
	const auto& i_oldCrvDigiMC = oldCrvDigiMCs.at(i_old_digi_mc);
	unsigned int i_new_digi_mc = i_old_digi_mc;
	if (_crvDigiMCIndexMapTag != "") {
	  if (_crvDigiMCIndexMap.checkInMap(i_old_digi_mc)) {
	    i_new_digi_mc = _crvDigiMCIndexMap.getCondensedIndex(i_old_digi_mc);
	  }
	  else {
	    continue; // to next old digi MC, since this one was compressed out...
	  }
	}
	const auto& i_newCrvDigiMC = newCrvDigiMCs.at(i_new_digi_mc);

	if (i_oldCrvDigiMC.GetStepPoints().size() > 0 && i_newCrvDigiMC.GetStepPoints().size()>0) {
	  const auto& i_oldStepPointMC = *i_oldCrvDigiMC.GetStepPoints().begin();
	  if (!i_oldStepPointMC.isAvailable()) {
	    continue; // this is a null step point
	  }
	  const auto& i_newStepPointMC = *i_newCrvDigiMC.GetStepPoints().begin();

	  const auto& i_old_digi_mc_barIndex = i_oldCrvDigiMC.GetScintillatorBarIndex();
	  const auto& i_new_digi_mc_barIndex = i_newCrvDigiMC.GetScintillatorBarIndex();
	  if (i_old_digi_mc_barIndex != i_new_digi_mc_barIndex) {
	    throw cet::exception("CompressDigiMCsCheck") << "Old and new CrvDigiMC's ScintillatorBarIndexs do not match" << std::endl;
	  }
	
	  const auto& i_old_step_barIndex = i_oldStepPointMC->barIndex();
	  const auto& i_new_step_barIndex = i_newStepPointMC->barIndex();
	  if (i_old_step_barIndex != i_new_step_barIndex) {
	    throw cet::exception("CompressDigiMCsCheck") << "Old and new CrvDigiMC's StepPointMC's BarIndexs do not match" << std::endl;
	  }
	
	  if (i_new_step_barIndex != i_new_digi_mc_barIndex) {
	    throw cet::exception("CompressDigiMCsCheck") << "New CrvDigiMC's BarIndex is inconsistent with its StepPointMC's BarIndex" << std::endl;
	  }
	
	  double old_time = _oldTOff.timeWithOffsetsApplied(*i_oldStepPointMC);
	  double new_time = _newTOff.timeWithOffsetsApplied(*i_newStepPointMC);
	  if (std::fabs(old_time - new_time) > 1e-5) {
	    throw cet::exception("CompressDigiMCsCheck") << "Old and new StepPointMC times with offsets applied do not match (CrvDigiMC)" << std::endl;
	  }
	}
      }
    }


    //////////////////////////////////////
    // Check CaloShowerSims
    art::Handle<CaloShowerSimCollection> oldCaloShowerSimHandle;
    event.getByLabel(_oldCaloShowerSimTag, oldCaloShowerSimHandle);

    art::Handle<CaloShowerSimCollection> newCaloShowerSimHandle;
    event.getByLabel(_newCaloShowerSimTag, newCaloShowerSimHandle);

    if (oldCaloShowerSimHandle.isValid() && newCaloShowerSimHandle.isValid()) {
      const auto& oldCaloShowerSims = *oldCaloShowerSimHandle;
      const auto& newCaloShowerSims = *newCaloShowerSimHandle;

      unsigned int n_old_calo_shower_sims = oldCaloShowerSims.size();
      for (unsigned int i_old_calo_shower_sim = 0; i_old_calo_shower_sim < n_old_calo_shower_sims; ++i_old_calo_shower_sim) {
	unsigned int i_new_calo_shower_sim = i_old_calo_shower_sim;
	const auto& i_oldCaloShowerSim = oldCaloShowerSims.at(i_old_calo_shower_sim);
	const auto& i_newCaloShowerSim = newCaloShowerSims.at(i_new_calo_shower_sim);

	const auto& i_oldCaloShowerStepPtr = *(i_oldCaloShowerSim.caloShowerSteps().begin());
	const auto& i_newCaloShowerStepPtr = *(i_newCaloShowerSim.caloShowerSteps().begin());
	if (i_oldCaloShowerStepPtr->momentumIn() != i_newCaloShowerStepPtr->momentumIn()) {
	  throw cet::exception("CompressDigiMCsCheck") << "Old and new CaloShowerStepPtrs do not match" << std::endl;
	}
      }
    }
  } // analyze(event)

  //================================================================


} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CompressDigiMCsCheck);
