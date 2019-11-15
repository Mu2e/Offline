////////////////////////////////////////////////////////////////////////
// Class:       StepPointsInDigis
// Plugin Type: producer (art v2_06_02)
// File:        StepPointsInDigis_module.cc
//
// Creates new StrawDigiMC and CrvDigiMC collections after creating new
// StepPointMC, SimParticle, GenParticle and SimParticleTimeMaps with all 
// unnecessary MC objects removed
//
// Generated at Wed Apr 12 16:10:46 2017 by Andrew Edmonds using cetskelgen
// from cetlib version v2_02_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include <memory>

#include "TTree.h"

#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/CrvDigiMCCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

#include "MCDataProducts/inc/GenId.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

namespace mu2e {
  class StepPointsInDigis;
}


class mu2e::StepPointsInDigis : public art::EDAnalyzer {
public:
  explicit StepPointsInDigis(fhicl::ParameterSet const & pset);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  StepPointsInDigis(StepPointsInDigis const &) = delete;
  StepPointsInDigis(StepPointsInDigis &&) = delete;
  StepPointsInDigis & operator = (StepPointsInDigis const &) = delete;
  StepPointsInDigis & operator = (StepPointsInDigis &&) = delete;

  // Required functions.
  void analyze(art::Event const& event) override;

  // Other functions
  void fillTree(const mu2e::StepPointMC& old_step);
  void fillTree(const mu2e::StrawGasStep& old_step);


private:

  // art tags for the input collections
  art::InputTag _strawDigiMCTag;
  art::InputTag _crvDigiMCTag;

  // handles to the old collections
  art::Handle<StrawDigiMCCollection> _strawDigiMCsHandle;
  art::Handle<CrvDigiMCCollection> _crvDigiMCsHandle;


  int _diagLevel;
  TTree* _steps;
  double _stepX, _stepY, _stepZ, _stepRawTime, _stepOffsettedTime, _stepEDep;
  unsigned _stepGenId;
  unsigned _stepProductId;

  TTree* _digis;
  double _digiTime;
  unsigned _digiProductId;

  SimParticleTimeOffset _toff;
};


mu2e::StepPointsInDigis::StepPointsInDigis(fhicl::ParameterSet const & pset)
  : art::EDAnalyzer(pset),
    _strawDigiMCTag(pset.get<art::InputTag>("strawDigiMCTag")),
    _crvDigiMCTag(pset.get<art::InputTag>("crvDigiMCTag")),
    _diagLevel(pset.get<int>("diagLevel", 0)),
    _toff(pset.get<fhicl::ParameterSet>("TimeOffsets"))
{
  // Call appropriate produces<>() functions here.
  art::ServiceHandle<art::TFileService> tfs;
  _steps = tfs->make<TTree>("_steps", "Steps in Digis");
  _steps->Branch("x", &_stepX);
  _steps->Branch("y", &_stepY);
  _steps->Branch("z", &_stepZ);
  _steps->Branch("raw_time", &_stepRawTime);
  _steps->Branch("offsetted_time", &_stepOffsettedTime);
  _steps->Branch("edep", &_stepEDep);
  _steps->Branch("genId", &_stepGenId);
  _steps->Branch("productId", &_stepProductId);

  _digis = tfs->make<TTree>("_digis", "Digis");
  _digis->Branch("time", &_digiTime);
  _digis->Branch("productId", &_digiProductId);
  _digis->Branch("step_raw_time", &_stepRawTime);
  _digis->Branch("step_offsetted_time", &_stepOffsettedTime);
}

void mu2e::StepPointsInDigis::analyze(art::Event const& event)
{
  // Implementation of required member function here.
  _toff.updateMap(event);

  event.getByLabel(_strawDigiMCTag, _strawDigiMCsHandle);
  const auto& strawDigiMCs = *_strawDigiMCsHandle;
  for (const auto& i_strawDigiMC : strawDigiMCs) {

    for(int i_end=0;i_end<StrawEnd::nends;++i_end){
      StrawEnd::End end = static_cast<StrawEnd::End>(i_end);
      auto const& old_step_point = i_strawDigiMC.strawGasStep(end);
      if (old_step_point.isAvailable()) {
	fillTree( *old_step_point );
      }
    }

    _digiProductId = _stepProductId;
    _digiTime = i_strawDigiMC.wireEndTime(StrawEnd::cal);
    _digis->Fill();
  }

  event.getByLabel(_crvDigiMCTag, _crvDigiMCsHandle);
  const auto& crvDigiMCs = *_crvDigiMCsHandle;
  for (const auto& i_crvDigiMC : crvDigiMCs) {

    for (const auto& i_step_mc : i_crvDigiMC.GetStepPoints()) {
      if (i_step_mc.isAvailable()) {
	fillTree(*i_step_mc);
      }
    }
  }
}

void mu2e::StepPointsInDigis::fillTree(const mu2e::StepPointMC& old_step) {

  _stepX = old_step.position().x();
  _stepY = old_step.position().y();
  _stepZ = old_step.position().z();
  _stepRawTime = old_step.time();
  _stepOffsettedTime = _toff.timeWithOffsetsApplied(old_step);
  _stepEDep = old_step.totalEDep();
  art::Ptr<SimParticle> simPtr = old_step.simParticle();
  while (simPtr->isSecondary()) {
    simPtr = simPtr->parent();
  }
  _stepGenId = simPtr->genParticle()->generatorId().id();
  _stepProductId = simPtr.id().value();
  _steps->Fill();
}

void mu2e::StepPointsInDigis::fillTree(const mu2e::StrawGasStep& old_step) {

  _stepX = old_step.position().x();
  _stepY = old_step.position().y();
  _stepZ = old_step.position().z();
  _stepRawTime = old_step.time();
  _stepOffsettedTime = _toff.timeWithOffsetsApplied(old_step);
  _stepEDep = old_step.totalEDep();
  art::Ptr<SimParticle> simPtr = old_step.simParticle();
  while (simPtr->isSecondary()) {
    simPtr = simPtr->parent();
  }
  _stepGenId = simPtr->genParticle()->generatorId().id();
  _stepProductId = simPtr.id().value();
  _steps->Fill();
}


DEFINE_ART_MODULE(mu2e::StepPointsInDigis)
