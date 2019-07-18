////////////////////////////////////////////////////////////////////////
// Class:       CompressStepPointMCs
// Plugin Type: producer (art v2_06_02)
// File:        CompressStepPointMCs_module.cc
//
// Takes an input StepPointMCCollection and removes StepPointMCs that
// don't pass an edep and time cut
//
// This module also removes any extraneous MC information (e.g. SimParticles,
// GenParticles) that are no longer pointed to
//
// A. Edmonds July 2017
//
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
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

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "Mu2eUtilities/inc/compressSimParticleCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"

#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"

namespace mu2e {
  class CompressStepPointMCs;

  class SimParticleSelector {
  public:
    SimParticleSelector() { }
    
    void push_back(cet::map_vector_key key) {
      m_keys.insert(key);
    }

    bool operator[]( cet::map_vector_key key ) const {
      return m_keys.find(key) != m_keys.end();
    }

    const std::set<cet::map_vector_key>& keys() const {
      return m_keys;
    }

    void clear() {
      m_keys.clear();
    }

  private:
    std::set<cet::map_vector_key> m_keys;
    
  };

  struct StepDiagInfo {
    Int_t _eventid;
    Float_t _stepRawTime;
    Float_t _stepTime;
    Float_t _stepEdep;
    Int_t _filtered;
    Int_t _nSimGenerations;

    static std::string leafNames() { return "eventid/I:raw_time/F:time/F:edep/F:filtered/I:nSimGenerations/I"; }
  };
}


class mu2e::CompressStepPointMCs : public art::EDFilter {
public:
  explicit CompressStepPointMCs(fhicl::ParameterSet const & pset);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CompressStepPointMCs(CompressStepPointMCs const &) = delete;
  CompressStepPointMCs(CompressStepPointMCs &&) = delete;
  CompressStepPointMCs & operator = (CompressStepPointMCs const &) = delete;
  CompressStepPointMCs & operator = (CompressStepPointMCs &&) = delete;

  // Required functions.
  bool filter(art::Event & event) override;

  // Other functions
  virtual void endJob() override;
  art::Ptr<StepPointMC> copyStepPointMC(const mu2e::StepPointMC& old_step);
  art::Ptr<CaloShowerStep> copyCaloShowerStep(const mu2e::CaloShowerStep& old_step);

private:

  // art tags for the input collections
  std::vector<art::InputTag> _stepPointMCTags;
  std::vector<art::InputTag> _caloShowerStepTags;
  art::InputTag _simParticleTag;
  std::vector<art::InputTag> _timeMapTags;

  // handles to the old collections
  art::Handle<StepPointMCCollection> _stepPointMCsHandle;
  art::Handle<CaloShowerStepCollection> _caloShowerStepsHandle;
  std::vector<SimParticleTimeMap> _oldTimeMaps;

  // unique_ptrs to the new output collections
  std::vector<std::unique_ptr<StepPointMCCollection> > _newStepPointMCs;
  std::vector<std::unique_ptr<CaloShowerStepCollection> > _newCaloShowerSteps;
  std::unique_ptr<SimParticleCollection> _newSimParticles;
  std::unique_ptr<GenParticleCollection> _newGenParticles;
  std::vector<std::unique_ptr<SimParticleTimeMap> > _newSimParticleTimeMaps;

  // for StepPointMCs, SimParticles and GenParticles we also need reference their new locations with art::Ptrs and so need their ProductIDs and Getters
  art::ProductID _newStepPointMCsPID;
  const art::EDProductGetter* _newStepPointMCGetter;
  art::ProductID _newCaloShowerStepsPID;
  const art::EDProductGetter* _newCaloShowerStepGetter;
  art::ProductID _newSimParticlesPID;
  const art::EDProductGetter* _oldSimParticleGetter;
  const art::EDProductGetter* _newSimParticleGetter;
  art::ProductID _newGenParticlesPID;
  const art::EDProductGetter* _newGenParticleGetter;

  // record the SimParticles that we are keeping so we can use compressSimParticleCollection to do all the work for us
  SimParticleSelector _simParticlesToKeep;

  SimParticleTimeOffset _toff;

  double _minTime;
  double _maxTime;
  double _minEdep;
  double _maxEdep;
  uint16_t _allStraw; // minimum straw # to keep all hits
  std::vector<uint16_t> _allPlanes; // planes in which to keep all hits

  int _diagLevel;
  TTree* _filterDiag;
  StepDiagInfo _stepDiag;
  TTree* _filterDiagCalo;
  StepDiagInfo _caloStepDiag;

  double _mbtime; // period of 1 microbunch

  unsigned _numInputEvents;
  unsigned _numOutputEvents;
};


mu2e::CompressStepPointMCs::CompressStepPointMCs(fhicl::ParameterSet const & pset)
  : art::EDFilter{pset},
    _stepPointMCTags(pset.get<std::vector<art::InputTag> >("stepPointMCTags")),
    _caloShowerStepTags(pset.get<std::vector<art::InputTag> >("caloShowerStepTags")),
    _simParticleTag(pset.get<art::InputTag>("simParticleTag")),
    _timeMapTags(pset.get<std::vector<art::InputTag> >("timeMapTags")),
    _toff(_timeMapTags),
    _minTime(pset.get<double>("minTime")),
    _maxTime(pset.get<double>("maxTime")),
    _minEdep(pset.get<double>("minEdep")),
    _maxEdep(pset.get<double>("maxEdep")),
    _allStraw(pset.get<uint16_t>("AllHitsStraw",90)), 
    _allPlanes(pset.get<std::vector<uint16_t>>("AllHitsPlanes",std::vector<uint16_t>{})), // planes to read all hits
    _diagLevel(pset.get<int>("diagLevel")),
    _numInputEvents(0), _numOutputEvents(0)
{
  // Call appropriate produces<>() functions here.
  for (std::vector<art::InputTag>::const_iterator i_tag = _stepPointMCTags.begin(); i_tag != _stepPointMCTags.end(); ++i_tag) {
    produces<StepPointMCCollection>( (*i_tag).instance() );
  }
  for (std::vector<art::InputTag>::const_iterator i_tag = _caloShowerStepTags.begin(); i_tag != _caloShowerStepTags.end(); ++i_tag) {
    produces<CaloShowerStepCollection>( (*i_tag).label() );
  }
  produces<SimParticleCollection>();
  produces<GenParticleCollection>();

  for (std::vector<art::InputTag>::const_iterator i_tag = _timeMapTags.begin(); i_tag != _timeMapTags.end(); ++i_tag) {
    produces<SimParticleTimeMap>( (*i_tag).label() );
  }

  if (_diagLevel > 0) {
    art::ServiceHandle<art::TFileService> tfs;
    _filterDiag = tfs->make<TTree>("_fdiag", "Filter Diag");
    _filterDiag->Branch("step", &_stepDiag, StepDiagInfo::leafNames().c_str());

    _filterDiagCalo = tfs->make<TTree>("_fdiagCalo", "Filter Diag");
    _filterDiagCalo->Branch("calo_step", &_caloStepDiag, StepDiagInfo::leafNames().c_str());
  }
}

bool mu2e::CompressStepPointMCs::filter(art::Event & event)
{
  // Implementation of required member function here.
  bool passed = false;
  _simParticlesToKeep.clear();
    
  _newSimParticles = std::unique_ptr<SimParticleCollection>(new SimParticleCollection);
  _newSimParticlesPID = event.getProductID<SimParticleCollection>();
  _newSimParticleGetter = event.productGetter(_newSimParticlesPID);
    
  _newGenParticles = std::unique_ptr<GenParticleCollection>(new GenParticleCollection);
  _newGenParticlesPID = event.getProductID<GenParticleCollection>();
  _newGenParticleGetter = event.productGetter(_newGenParticlesPID);

  _oldTimeMaps.clear();
  for (std::vector<art::InputTag>::const_iterator i_tag = _timeMapTags.begin(); i_tag != _timeMapTags.end(); ++i_tag) {
    art::Handle<SimParticleTimeMap> i_timeMapHandle;
    event.getByLabel(*i_tag, i_timeMapHandle);

    if (!i_timeMapHandle.isValid()) {
      throw cet::exception("CompressStepPointMCs") << "Couldn't find SimParticleTimeMap " << *i_tag << " in event\n";
    }
    _oldTimeMaps.push_back(*i_timeMapHandle);
  }

  _newSimParticleTimeMaps.clear();
  for (std::vector<art::InputTag>::const_iterator i_tag = _timeMapTags.begin(); i_tag != _timeMapTags.end(); ++i_tag) {
    _newSimParticleTimeMaps.push_back(std::unique_ptr<SimParticleTimeMap>(new SimParticleTimeMap));
  }

  ConditionsHandle<AcceleratorParams> accPar("ignored");
  _mbtime = accPar->deBuncherPeriod;
  _toff.updateMap(event);

  _newStepPointMCs.clear();
  for (const auto& i_stepTag : _stepPointMCTags) {
    _newStepPointMCs.push_back(std::unique_ptr<StepPointMCCollection>(new StepPointMCCollection));  
    _newStepPointMCsPID = event.getProductID<StepPointMCCollection>( i_stepTag.instance() );
    _newStepPointMCGetter = event.productGetter(_newStepPointMCsPID);

    event.getByLabel(i_stepTag, _stepPointMCsHandle);
    const auto& stepPointMCs = *_stepPointMCsHandle;
    for (const auto& i_stepPointMC : stepPointMCs) {
      double i_edep = i_stepPointMC.totalEDep();
      double i_time = std::fmod(_toff.timeWithOffsetsApplied(i_stepPointMC), _mbtime);
      while (i_time < 0) {
	i_time += _mbtime;
      }
      
      if (_diagLevel > 0) {
	_stepDiag._eventid = event.id().event();
	_stepDiag._stepRawTime = i_stepPointMC.time();
	_stepDiag._stepTime = i_time;
	_stepDiag._stepEdep = i_edep;
	_stepDiag._filtered = 0;
      }
     
     // check if we are keeping all hits for this straw
      StrawId sid = i_stepPointMC.strawId();
      bool keepall = sid.straw() >= _allStraw &&
	(std::find(_allPlanes.begin(),_allPlanes.end(),sid.plane()) != _allPlanes.end());

      if ( keepall || 
	  ((i_edep > _minEdep && i_edep < _maxEdep) &&
	   (i_time > _minTime && i_time < _maxTime) ) )  {

	if (_diagLevel > 0) {
	  _stepDiag._filtered = 1;
	}
	if(_diagLevel > 1 && keepall)
	  std::cout << " Keeping hit in straw " << sid << " time " << i_time << std::endl;
	copyStepPointMC(i_stepPointMC);
      }
      
      if (_diagLevel > 0) {
	_filterDiag->Fill();
      }
    }
  }

  _newCaloShowerSteps.clear();
  for (const auto& i_caloStepTag : _caloShowerStepTags) {
    _newCaloShowerSteps.push_back(std::unique_ptr<CaloShowerStepCollection>(new CaloShowerStepCollection));  
    _newCaloShowerStepsPID = event.getProductID<CaloShowerStepCollection>( i_caloStepTag.label() );
    _newCaloShowerStepGetter = event.productGetter(_newCaloShowerStepsPID);

    event.getByLabel(i_caloStepTag, _caloShowerStepsHandle);
    const auto& caloShowerSteps = *_caloShowerStepsHandle;
    for (const auto& i_caloShowerStep : caloShowerSteps) {

      double i_edep = i_caloShowerStep.energyMC();
      double i_time;
      if (i_caloShowerStep.simParticle().get()) { // see note below in copyCaloShowerStep
	i_time = std::fmod(i_caloShowerStep.timeStepMC()+_toff.totalTimeOffset(i_caloShowerStep.simParticle()), _mbtime);
      }
      else {
	continue;
      }
      
      while (i_time < 0) {
	i_time += _mbtime;
      }


      if (_diagLevel > 0) {
	_caloStepDiag._eventid = event.id().event();
	_caloStepDiag._stepRawTime = i_caloShowerStep.timeStepMC();
	_caloStepDiag._stepTime = i_time;
	_caloStepDiag._stepEdep = i_edep;
	_caloStepDiag._filtered = 0;
      }

      if ( (i_edep > _minEdep && i_edep < _maxEdep) &&
	   (i_time > _minTime && i_time < _maxTime) ) {
	
	if (_diagLevel > 0) {
	  _caloStepDiag._filtered = 1;
	}
	copyCaloShowerStep(i_caloShowerStep);
      }

      if (_diagLevel > 0) {
	_filterDiagCalo->Fill();
      }

    }
  }

  
  // Now compress the SimParticleCollections into their new collections
  const auto& oldSimParticles = event.getValidHandle<SimParticleCollection>(_simParticleTag);
  art::ProductID i_product_id = oldSimParticles.id();
  compressSimParticleCollection(_newSimParticlesPID, _newSimParticleGetter, *oldSimParticles, 
				_simParticlesToKeep, *_newSimParticles);

  passed = !_newSimParticles->empty();
  for(auto& i : *_newSimParticles) {
      
    mu2e::SimParticle& newsim = i.second;
    if(!newsim.genParticle().isNull()) { // will crash if not resolvable

      // Copy GenParticle to the new collection
      _newGenParticles->emplace_back(*newsim.genParticle());
      newsim.genParticle() = art::Ptr<GenParticle>(_newGenParticlesPID, _newGenParticles->size()-1, _newGenParticleGetter);
    }
      
    // Update the time maps
    art::Ptr<SimParticle> oldSimPtr(i_product_id, newsim.id().asUint(), _oldSimParticleGetter);
    art::Ptr<SimParticle> newSimPtr(_newSimParticlesPID, newsim.id().asUint(), _newSimParticleGetter);
    for (std::vector<SimParticleTimeMap>::const_iterator i_time_map = _oldTimeMaps.begin(); i_time_map != _oldTimeMaps.end(); ++i_time_map) {
      size_t i_element = i_time_map - _oldTimeMaps.begin();
      
      const SimParticleTimeMap& i_oldTimeMap = *i_time_map;
      SimParticleTimeMap& i_newTimeMap = *_newSimParticleTimeMaps.at(i_element);
      auto it = i_oldTimeMap.find(oldSimPtr);
      if (it != i_oldTimeMap.end()) {
	i_newTimeMap[newSimPtr] = it->second;
      }
    }
  }

  // Now add everything to the event
  for (std::vector<art::InputTag>::const_iterator i_tag = _stepPointMCTags.begin(); i_tag != _stepPointMCTags.end(); ++i_tag) {
    size_t i_element = i_tag - _stepPointMCTags.begin();
    event.put(std::move(_newStepPointMCs.at(i_element)), (*i_tag).instance());
  }
  for (std::vector<art::InputTag>::const_iterator i_tag = _caloShowerStepTags.begin(); i_tag != _caloShowerStepTags.end(); ++i_tag) {
    size_t i_element = i_tag - _caloShowerStepTags.begin();
    event.put(std::move(_newCaloShowerSteps.at(i_element)), (*i_tag).label());
  }
  event.put(std::move(_newSimParticles));
  event.put(std::move(_newGenParticles));

  for (std::vector<art::InputTag>::const_iterator i_tag = _timeMapTags.begin(); i_tag != _timeMapTags.end(); ++i_tag) {
    size_t i_element = i_tag - _timeMapTags.begin();
    event.put(std::move(_newSimParticleTimeMaps.at(i_element)), (*i_tag).label());
  }

  ++_numInputEvents;
  if (passed) { ++_numOutputEvents; }
  return passed;
}

art::Ptr<mu2e::StepPointMC> mu2e::CompressStepPointMCs::copyStepPointMC(const mu2e::StepPointMC& old_step) {

  _simParticlesToKeep.push_back(old_step.simParticle()->id());
  art::Ptr<SimParticle> newSimPtr(_newSimParticlesPID, old_step.simParticle()->id().asUint(), _newSimParticleGetter);

  // Also need to add all the parents (and change their genParticles) too
  art::Ptr<SimParticle> childPtr = old_step.simParticle();
  art::Ptr<SimParticle> parentPtr = childPtr->parent();

  _stepDiag._nSimGenerations = 1;
  while (parentPtr) {
    _simParticlesToKeep.push_back(parentPtr->id());

    childPtr = parentPtr;
    parentPtr = parentPtr->parent();
    ++_stepDiag._nSimGenerations;
  }
  
  StepPointMC new_step(old_step);
  new_step.simParticle() = newSimPtr;

  _newStepPointMCs.back()->push_back(new_step);
  
  return art::Ptr<StepPointMC>(_newStepPointMCsPID, _newStepPointMCs.back()->size()-1, _newStepPointMCGetter);
}

art::Ptr<mu2e::CaloShowerStep> mu2e::CompressStepPointMCs::copyCaloShowerStep(const mu2e::CaloShowerStep& old_step) {

  // Need this if-statement because sometimes the SimParticle that is being Ptr'd to 
  // is not there... The Ptr itself is valid (i.e. old_step.simParticle().isNonnull() returns true)
  // but there is no object there and so when we try to get the id of the SimParticle
  // there is a segfault
  if (old_step.simParticle().get()) {
    _simParticlesToKeep.push_back(old_step.simParticle()->id());
    art::Ptr<SimParticle> newSimPtr(_newSimParticlesPID, old_step.simParticle()->id().asUint(), _newSimParticleGetter);

    // Also need to add all the parents (and change their genParticles) too
    art::Ptr<SimParticle> childPtr = old_step.simParticle();
    art::Ptr<SimParticle> parentPtr = childPtr->parent();

    _caloStepDiag._nSimGenerations = 1;  
    while (parentPtr) {
      _simParticlesToKeep.push_back(parentPtr->id());
      
      childPtr = parentPtr;
      parentPtr = parentPtr->parent();
      ++_caloStepDiag._nSimGenerations;
    }
    
    CaloShowerStep new_step(old_step.volumeId(), newSimPtr, old_step.nCompress(), old_step.timeStepMC(), old_step.energyMC(), old_step.momentumIn(), old_step.positionIn(), old_step.position(), old_step.covPosition());
    
    _newCaloShowerSteps.back()->push_back(new_step);
    
    return art::Ptr<CaloShowerStep>(_newCaloShowerStepsPID, _newCaloShowerSteps.back()->size()-1, _newCaloShowerStepGetter);
  }
  else {
    return art::Ptr<CaloShowerStep>();
  }
}

void mu2e::CompressStepPointMCs::endJob() {

  mf::LogInfo("Summary")
    << "CompressStepPointMCs_module stats: passed = " << _numOutputEvents << " / " << _numInputEvents;
}


DEFINE_ART_MODULE(mu2e::CompressStepPointMCs)
