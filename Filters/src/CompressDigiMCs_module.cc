////////////////////////////////////////////////////////////////////////
// Class:       CompressDigiMCs
// Plugin Type: producer (art v2_06_02)
// File:        CompressDigiMCs_module.cc
//
// Creates new StrawDigiMC and CrvDigiMC collections after creating new
// StepPointMC, SimParticle, GenParticle and SimParticleTimeMaps with all 
// unnecessary MC objects removed.
//
// Also creates new CaloShowerStep, CaloShowerStepRO and CaloShowerSim collections after
// remapping the art::Ptrs to the SimParticles
//
// Generated at Wed Apr 12 16:10:46 2017 by Andrew Edmonds using cetskelgen
// from cetlib version v2_02_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include <memory>

#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/CrvDigiMCCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/CaloShowerSimCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepROCollection.hh"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "Mu2eUtilities/inc/compressSimParticleCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"
#include "MCDataProducts/inc/SimParticleRemapping.hh"

namespace mu2e {
  class CompressDigiMCs;

  typedef std::set<art::Ptr<SimParticle> > SimParticleSet;

  class SimParticleSelector {
  public:
    SimParticleSelector(const SimParticleSet& simPartSet) { 
      for (const auto& i_simPart : simPartSet) {
	cet::map_vector_key key = cet::map_vector_key(i_simPart.key());
	m_keys.insert(key);
      }
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

  typedef std::map<art::Ptr<mu2e::CaloShowerStep>, art::Ptr<mu2e::CaloShowerStep> > CaloShowerStepRemap;
  typedef std::string InstanceLabel;
  typedef std::map<cet::map_vector_key, cet::map_vector_key> KeyRemap;
}


class mu2e::CompressDigiMCs : public art::EDProducer {
public:
  explicit CompressDigiMCs(fhicl::ParameterSet const & pset);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CompressDigiMCs(CompressDigiMCs const &) = delete;
  CompressDigiMCs(CompressDigiMCs &&) = delete;
  CompressDigiMCs & operator = (CompressDigiMCs const &) = delete;
  CompressDigiMCs & operator = (CompressDigiMCs &&) = delete;

  // Required functions.
  void produce(art::Event & event) override;

  // Other functions
  void copyStrawDigiMC(const mu2e::StrawDigiMC& old_straw_digi_mc);
  void copyCrvDigiMC(const mu2e::CrvDigiMC& old_crv_digi_mc);
  art::Ptr<StepPointMC> copyStepPointMC(const mu2e::StepPointMC& old_step, const InstanceLabel& instance);
  art::Ptr<mu2e::CaloShowerStep> copyCaloShowerStep(const mu2e::CaloShowerStep& old_calo_shower_step);
  void copyCaloShowerSim(const mu2e::CaloShowerSim& old_calo_shower_sim, const CaloShowerStepRemap& remap);
  void copyCaloShowerStepRO(const mu2e::CaloShowerStepRO& old_calo_shower_step_ro, const CaloShowerStepRemap& remap);
  void keepSimParticle(const art::Ptr<SimParticle>& sim_ptr);

private:

  // art tags for the input collections
  art::InputTag _strawDigiMCTag;
  art::InputTag _crvDigiMCTag;

  std::vector<art::InputTag> _simParticleTags;
  std::vector<art::InputTag> _extraStepPointMCTags;
  std::vector<art::InputTag> _timeMapTags;
  std::vector<art::InputTag> _caloShowerStepTags;
  art::InputTag _caloShowerSimTag;
  art::InputTag _caloShowerStepROTag;

  // handles to the old collections
  art::Handle<StrawDigiMCCollection> _strawDigiMCsHandle;
  art::Handle<CrvDigiMCCollection> _crvDigiMCsHandle;
  std::vector<SimParticleTimeMap> _oldTimeMaps;
  art::Handle<CaloShowerStepCollection> _caloShowerStepsHandle;
  art::Handle<CaloShowerSimCollection> _caloShowerSimsHandle;
  art::Handle<CaloShowerStepROCollection> _caloShowerStepROsHandle;

  // unique_ptrs to the new output collections
  std::unique_ptr<StrawDigiMCCollection> _newStrawDigiMCs;
  std::unique_ptr<CrvDigiMCCollection> _newCrvDigiMCs;
  std::map<InstanceLabel, std::unique_ptr<StepPointMCCollection> > _newStepPointMCs;
  std::unique_ptr<SimParticleCollection> _newSimParticles;
  std::unique_ptr<GenParticleCollection> _newGenParticles;
  std::vector<std::unique_ptr<SimParticleTimeMap> > _newSimParticleTimeMaps;
  std::unique_ptr<CaloShowerStepCollection> _newCaloShowerSteps;
  std::unique_ptr<CaloShowerSimCollection> _newCaloShowerSims;
  std::unique_ptr<CaloShowerStepROCollection> _newCaloShowerStepROs;

  // for StepPointMCs, SimParticles and GenParticles we also need reference their new locations with art::Ptrs and so need their ProductIDs and Getters
  std::map<InstanceLabel, art::ProductID> _newStepPointMCsPID;
  std::map<InstanceLabel, const art::EDProductGetter*> _newStepPointMCGetter;
  art::ProductID _newSimParticlesPID;
  const art::EDProductGetter* _newSimParticleGetter;
  art::ProductID _newGenParticlesPID;
  const art::EDProductGetter* _newGenParticleGetter;
  art::ProductID _newCaloShowerStepsPID;
  const art::EDProductGetter* _newCaloShowerStepGetter;
  std::map<art::ProductID, const art::EDProductGetter*> _oldCaloShowerStepGetter;

  // record the SimParticles that we are keeping so we can use compressSimParticleCollection to do all the work for us
  std::map<art::ProductID, SimParticleSet> _simParticlesToKeep;

  InstanceLabel _trackerOutputInstanceLabel;
  InstanceLabel _crvOutputInstanceLabel;
  std::vector<InstanceLabel> _newStepPointMCInstances;
};


mu2e::CompressDigiMCs::CompressDigiMCs(fhicl::ParameterSet const & pset)
  : _strawDigiMCTag(pset.get<art::InputTag>("strawDigiMCTag")),
    _crvDigiMCTag(pset.get<art::InputTag>("crvDigiMCTag")),
    _simParticleTags(pset.get<std::vector<art::InputTag> >("simParticleTags")),
    _extraStepPointMCTags(pset.get<std::vector<art::InputTag> >("extraStepPointMCTags")),
    _timeMapTags(pset.get<std::vector<art::InputTag> >("timeMapTags")),
    _caloShowerStepTags(pset.get<std::vector<art::InputTag> >("caloShowerStepTags")),
    _caloShowerSimTag(pset.get<art::InputTag>("caloShowerSimTag")),
    _caloShowerStepROTag(pset.get<art::InputTag>("caloShowerStepROTag")),
    _trackerOutputInstanceLabel(pset.get<std::string>("trackerOutputInstanceLabel", "tracker")),
    _crvOutputInstanceLabel(pset.get<std::string>("crvOutputInstanceLabel", "CRV"))
{
  // Call appropriate produces<>() functions here.
  produces<StrawDigiMCCollection>();
  produces<CrvDigiMCCollection>();

  _newStepPointMCInstances.push_back(_trackerOutputInstanceLabel); // will always be a tracker and CRV instance
  _newStepPointMCInstances.push_back("CRV");     // filled with the StepPointMCs referenced by their DigiMCs
  for (std::vector<art::InputTag>::const_iterator i_tag = _extraStepPointMCTags.begin(); i_tag != _extraStepPointMCTags.end(); ++i_tag) {
    _newStepPointMCInstances.push_back( (*i_tag).instance() );
  }

  for (const auto& i_instance : _newStepPointMCInstances) {
    produces<StepPointMCCollection>( i_instance );
  }

  produces<SimParticleCollection>();
  produces<GenParticleCollection>();

  for (std::vector<art::InputTag>::const_iterator i_tag = _timeMapTags.begin(); i_tag != _timeMapTags.end(); ++i_tag) {
    produces<SimParticleTimeMap>( (*i_tag).label() );
  }

  produces<CaloShowerStepCollection>();
  produces<CaloShowerSimCollection>();
  produces<CaloShowerStepROCollection>();
}

void mu2e::CompressDigiMCs::produce(art::Event & event)
{
  // Implementation of required member function here.
  _newStrawDigiMCs = std::unique_ptr<StrawDigiMCCollection>(new StrawDigiMCCollection);  
  _newCrvDigiMCs = std::unique_ptr<CrvDigiMCCollection>(new CrvDigiMCCollection);  

  for (const auto& i_instance : _newStepPointMCInstances) {
    _newStepPointMCs[i_instance] = std::unique_ptr<StepPointMCCollection>(new StepPointMCCollection);
    _newStepPointMCsPID[i_instance] = getProductID<StepPointMCCollection>(i_instance);
    _newStepPointMCGetter[i_instance] = event.productGetter(_newStepPointMCsPID[i_instance]);
  }

  _newSimParticles = std::unique_ptr<SimParticleCollection>(new SimParticleCollection);
  _newSimParticlesPID = getProductID<SimParticleCollection>();
  _newSimParticleGetter = event.productGetter(_newSimParticlesPID);
    
  _newGenParticles = std::unique_ptr<GenParticleCollection>(new GenParticleCollection);
  _newGenParticlesPID = getProductID<GenParticleCollection>();
  _newGenParticleGetter = event.productGetter(_newGenParticlesPID);

  // Create all the new collections, ProductIDs and product getters for the SimParticles and GenParticles
  // There is one for each background frame plus one for the primary event
  for (std::vector<art::InputTag>::const_iterator i_tag = _simParticleTags.begin(); i_tag != _simParticleTags.end(); ++i_tag) {
    const auto& oldSimParticles = event.getValidHandle<SimParticleCollection>(*i_tag);
    art::ProductID i_product_id = oldSimParticles.id();
    
    _simParticlesToKeep[i_product_id].clear();
  }

  _oldTimeMaps.clear();
  for (std::vector<art::InputTag>::const_iterator i_tag = _timeMapTags.begin(); i_tag != _timeMapTags.end(); ++i_tag) {
    art::Handle<SimParticleTimeMap> i_timeMapHandle;
    event.getByLabel(*i_tag, i_timeMapHandle);

    if (!i_timeMapHandle.isValid()) {
      throw cet::exception("CompressDigiMCs") << "Couldn't find SimParticleTimeMap " << *i_tag << " in event\n";
    }
    _oldTimeMaps.push_back(*i_timeMapHandle);
  }

  _newSimParticleTimeMaps.clear();
  for (std::vector<art::InputTag>::const_iterator i_tag = _timeMapTags.begin(); i_tag != _timeMapTags.end(); ++i_tag) {
    _newSimParticleTimeMaps.push_back(std::unique_ptr<SimParticleTimeMap>(new SimParticleTimeMap));
  }

  event.getByLabel(_strawDigiMCTag, _strawDigiMCsHandle);
  const auto& strawDigiMCs = *_strawDigiMCsHandle;
  for (const auto& i_strawDigiMC : strawDigiMCs) {
    copyStrawDigiMC(i_strawDigiMC);
  }

  if (_crvDigiMCTag != "") {
    event.getByLabel(_crvDigiMCTag, _crvDigiMCsHandle);
    const auto& crvDigiMCs = *_crvDigiMCsHandle;
    for (const auto& i_crvDigiMC : crvDigiMCs) {
      copyCrvDigiMC(i_crvDigiMC);
    }
  }

  CaloShowerStepRemap caloShowerStepRemap;
  _newCaloShowerSteps = std::unique_ptr<CaloShowerStepCollection>(new CaloShowerStepCollection);
  _newCaloShowerStepsPID = getProductID<CaloShowerStepCollection>();
  _newCaloShowerStepGetter = event.productGetter(_newCaloShowerStepsPID);
  for (std::vector<art::InputTag>::const_iterator i_tag = _caloShowerStepTags.begin(); i_tag != _caloShowerStepTags.end(); ++i_tag) {
    const auto& oldCaloShowerSteps = event.getValidHandle<CaloShowerStepCollection>(*i_tag);
    art::ProductID i_product_id = oldCaloShowerSteps.id();
    _oldCaloShowerStepGetter[i_product_id] = event.productGetter(i_product_id);

    for (CaloShowerStepCollection::const_iterator i_caloShowerStep = oldCaloShowerSteps->begin(); i_caloShowerStep != oldCaloShowerSteps->end(); ++i_caloShowerStep) {
      art::Ptr<mu2e::CaloShowerStep> oldShowerStepPtr(i_product_id,  i_caloShowerStep - oldCaloShowerSteps->begin(), _oldCaloShowerStepGetter[i_product_id]);
      art::Ptr<mu2e::CaloShowerStep> newShowerStepPtr = copyCaloShowerStep(*i_caloShowerStep);
      caloShowerStepRemap[oldShowerStepPtr] = newShowerStepPtr;
    }
  }

  _newCaloShowerSims = std::unique_ptr<CaloShowerSimCollection>(new CaloShowerSimCollection);
  event.getByLabel(_caloShowerSimTag, _caloShowerSimsHandle);
  const auto& caloShowerSims = *_caloShowerSimsHandle;
  for (const auto& i_caloShowerSim : caloShowerSims) {
    copyCaloShowerSim(i_caloShowerSim, caloShowerStepRemap);
  }

  _newCaloShowerStepROs = std::unique_ptr<CaloShowerStepROCollection>(new CaloShowerStepROCollection);
  event.getByLabel(_caloShowerStepROTag, _caloShowerStepROsHandle);
  const auto& caloShowerStepROs = *_caloShowerStepROsHandle;
  for (const auto& i_caloShowerStepRO : caloShowerStepROs) {
    copyCaloShowerStepRO(i_caloShowerStepRO, caloShowerStepRemap);
  }
  
  // Get the hits from the virtualdetector
  for (std::vector<art::InputTag>::const_iterator i_tag = _extraStepPointMCTags.begin(); i_tag != _extraStepPointMCTags.end(); ++i_tag) {
    const auto& stepPointMCs = event.getValidHandle<StepPointMCCollection>(*i_tag);
    for (const auto& stepPointMC : *stepPointMCs) {
      for (const auto& simPartsToKeep : _simParticlesToKeep) {
	const art::ProductID& oldProdID = simPartsToKeep.first;
	if (stepPointMC.simParticle().id() != oldProdID) {
	  continue;
	}
	const SimParticleSet& alreadyKeptSimParts = simPartsToKeep.second;
	for (const auto& alreadyKeptSimPart : alreadyKeptSimParts) {
	  if (stepPointMC.simParticle() == alreadyKeptSimPart) {
	    copyStepPointMC(stepPointMC, (*i_tag).instance() );
	  }
	}
      }
    }
  }
  
  // Now compress the SimParticleCollections into their new collections
  KeyRemap* keyRemap = new KeyRemap;
  SimParticleRemapping remap;
  unsigned int keep_size = 0;
  for (std::vector<art::InputTag>::const_iterator i_tag = _simParticleTags.begin(); i_tag != _simParticleTags.end(); ++i_tag) {
    keyRemap->clear();
    const auto& oldSimParticles = event.getValidHandle<SimParticleCollection>(*i_tag);
    art::ProductID i_product_id = oldSimParticles.id();
    SimParticleSelector simPartSelector(_simParticlesToKeep[i_product_id]);
    keep_size += _simParticlesToKeep[i_product_id].size();
    compressSimParticleCollection(_newSimParticlesPID, _newSimParticleGetter, *oldSimParticles, 
				  simPartSelector, *_newSimParticles, keyRemap);

    // Fill out the SimParticleRemapping
    for (const auto& i_keptSimPart : _simParticlesToKeep[i_product_id]) {
      cet::map_vector_key oldKey = cet::map_vector_key(i_keptSimPart.key());
      cet::map_vector_key newKey = keyRemap->at(oldKey);
      remap[i_keptSimPart] = art::Ptr<SimParticle>(_newSimParticlesPID, newKey.asUint(), _newSimParticleGetter);
    }
  }
  if (keep_size != _newSimParticles->size()) {
    throw cet::exception("CompressDigiMCs") << "Number of SimParticles in output collection (" << _newSimParticles->size() << ") does not match the number of SimParticles we wanted to keep (" << keep_size << ")" << std::endl;
  }

  // Loop through the new SimParticles to keep any GenParticles
  for (auto& i_simParticle : *_newSimParticles) {
    mu2e::SimParticle& newsim = i_simParticle.second;
    if(newsim.genParticle().isNonnull()) { // will crash if not resolvable
      
      // Copy GenParticle to the new collection
      _newGenParticles->emplace_back(*newsim.genParticle());
      newsim.genParticle() = art::Ptr<GenParticle>(_newGenParticlesPID, _newGenParticles->size()-1, _newGenParticleGetter);
    }
  }

  
  // Now update all objects with SimParticlePtrs
  // Update the time maps
  for (std::vector<SimParticleTimeMap>::const_iterator i_time_map = _oldTimeMaps.begin(); i_time_map != _oldTimeMaps.end(); ++i_time_map) {
    size_t i_element = i_time_map - _oldTimeMaps.begin();
	
    const SimParticleTimeMap& i_oldTimeMap = *i_time_map;
    SimParticleTimeMap& i_newTimeMap = *_newSimParticleTimeMaps.at(i_element);
    for (const auto& timeMapPair : i_oldTimeMap) {
      art::Ptr<SimParticle> oldSimPtr = timeMapPair.first;
      const auto& newSimPtrIter = remap.find(oldSimPtr);
      if (newSimPtrIter != remap.end()) {
	art::Ptr<SimParticle> newSimPtr = newSimPtrIter->second;
	i_newTimeMap[newSimPtr] = timeMapPair.second;
      }
    }
  }

  // Update the StepPointMCs
  for (const auto& i_instance : _newStepPointMCInstances) {
    for (auto& i_stepPointMC : *_newStepPointMCs.at(i_instance)) {
      art::Ptr<SimParticle> newSimPtr = remap.at(i_stepPointMC.simParticle());
      i_stepPointMC.simParticle() = newSimPtr;
    }
  }

  // Update the CaloShowerSteps
  for (auto& i_caloShowerStep : *_newCaloShowerSteps) {
    art::Ptr<SimParticle> newSimPtr = remap.at(i_caloShowerStep.simParticle());
    i_caloShowerStep.setSimParticle(newSimPtr);
  }

  // Update the CaloShowerSims
  for (auto& i_caloShowerSim : *_newCaloShowerSims) {
    art::Ptr<SimParticle> newSimPtr = remap.at(i_caloShowerSim.sim());
    i_caloShowerSim.setSimParticle(newSimPtr);
  }

  // Update the CrvDigiMCs
  for (auto& i_crvDigiMC : *_newCrvDigiMCs) {
    art::Ptr<SimParticle> oldSimPtr = i_crvDigiMC.GetSimParticle();
    art::Ptr<SimParticle> newSimPtr;
    if (oldSimPtr.isNonnull()) { // if the old CrvDigiMC doesn't have a null ptr for the SimParticle...
      newSimPtr = remap.at(oldSimPtr);
    }
    else {
      newSimPtr = art::Ptr<SimParticle>();
    }
    i_crvDigiMC.setSimParticle(newSimPtr);
  }


  // Now add everything to the event
  for (const auto& i_instance : _newStepPointMCInstances) {
    event.put(std::move(_newStepPointMCs.at(i_instance)), i_instance);
  }
  event.put(std::move(_newStrawDigiMCs));
  event.put(std::move(_newCrvDigiMCs));

  event.put(std::move(_newSimParticles));
  event.put(std::move(_newGenParticles));

  for (std::vector<art::InputTag>::const_iterator i_tag = _timeMapTags.begin(); i_tag != _timeMapTags.end(); ++i_tag) {
    size_t i_element = i_tag - _timeMapTags.begin();
    event.put(std::move(_newSimParticleTimeMaps.at(i_element)), (*i_tag).label());
  }

  event.put(std::move(_newCaloShowerSteps));
  event.put(std::move(_newCaloShowerSims));
  event.put(std::move(_newCaloShowerStepROs));
}

void mu2e::CompressDigiMCs::copyStrawDigiMC(const mu2e::StrawDigiMC& old_straw_digi_mc) {

  // Need to update the Ptrs for the StepPointMCs
  art::Ptr<StepPointMC> newTriggerStepPtr[StrawEnd::nends];
  for(int i_end=0;i_end<StrawEnd::nends;++i_end){
    StrawEnd::End end = static_cast<StrawEnd::End>(i_end);

    const art::Ptr<StepPointMC>& old_step_point = old_straw_digi_mc.stepPointMC(end);
    if (old_step_point.isAvailable()) {
      newTriggerStepPtr[i_end] = copyStepPointMC( *old_step_point, _trackerOutputInstanceLabel );
    }
  }
  
  std::vector<art::Ptr<StepPointMC> > newWaveformStepPtrs;
  for (const auto& i_step_mc : old_straw_digi_mc.stepPointMCs()) {
    if (i_step_mc.isAvailable()) {
      newWaveformStepPtrs.push_back(copyStepPointMC(*i_step_mc, _trackerOutputInstanceLabel));
    }
  }
  
  StrawDigiMC new_straw_digi_mc(old_straw_digi_mc, newTriggerStepPtr, newWaveformStepPtrs); // copy everything except the Ptrs from the old StrawDigiMC  
  _newStrawDigiMCs->push_back(new_straw_digi_mc);
}

void mu2e::CompressDigiMCs::copyCrvDigiMC(const mu2e::CrvDigiMC& old_crv_digi_mc) {

  // Need to update the Ptrs for the StepPointMCs
  std::vector<art::Ptr<StepPointMC> > newStepPtrs;
  for (const auto& i_step_mc : old_crv_digi_mc.GetStepPoints()) {
    if (i_step_mc.isAvailable()) {
      newStepPtrs.push_back(copyStepPointMC(*i_step_mc, _crvOutputInstanceLabel));
    }
  }

  CrvDigiMC new_crv_digi_mc(old_crv_digi_mc);
  new_crv_digi_mc.setStepPoints(newStepPtrs);

  _newCrvDigiMCs->push_back(new_crv_digi_mc);
}

art::Ptr<mu2e::CaloShowerStep> mu2e::CompressDigiMCs::copyCaloShowerStep(const mu2e::CaloShowerStep& old_calo_shower_step) {

  // Need this if-statement because sometimes the SimParticle that is being Ptr'd to 
  // is not there... The Ptr itself is valid (i.e. old_step.simParticle().isNonnull() returns true)
  // but there is no object there and so when we try to get the id of the SimParticle
  // there is a segfault
  if (old_calo_shower_step.simParticle().get()) {
    art::Ptr<SimParticle> oldSimPtr = old_calo_shower_step.simParticle();

    keepSimParticle(oldSimPtr);

    CaloShowerStep new_calo_shower_step = old_calo_shower_step;

    _newCaloShowerSteps->push_back(new_calo_shower_step);

    return art::Ptr<mu2e::CaloShowerStep>(_newCaloShowerStepsPID, _newCaloShowerSteps->size(), _newCaloShowerStepGetter);
  }
  else {
    return art::Ptr<CaloShowerStep>();
  }
}

void mu2e::CompressDigiMCs::copyCaloShowerSim(const mu2e::CaloShowerSim& old_calo_shower_sim, const CaloShowerStepRemap& remap) {

  art::Ptr<SimParticle> oldSimPtr = old_calo_shower_sim.sim();
  keepSimParticle(oldSimPtr);

  const auto& caloShowerStepPtrs = old_calo_shower_sim.caloShowerSteps();
  std::vector<art::Ptr<CaloShowerStep> > newCaloShowerStepPtrs;
  for (const auto& i_caloShowerStepPtr : caloShowerStepPtrs) {
    newCaloShowerStepPtrs.push_back(remap.at(i_caloShowerStepPtr));
  }

  CaloShowerSim new_calo_shower_sim = old_calo_shower_sim;
  new_calo_shower_sim.setCaloShowerSteps(newCaloShowerStepPtrs);

  _newCaloShowerSims->push_back(new_calo_shower_sim);
}

void mu2e::CompressDigiMCs::copyCaloShowerStepRO(const mu2e::CaloShowerStepRO& old_calo_shower_step_ro, const CaloShowerStepRemap& remap) {

  const auto& caloShowerStepPtr = old_calo_shower_step_ro.caloShowerStep();
  CaloShowerStepRO new_calo_shower_step_ro = old_calo_shower_step_ro;
  new_calo_shower_step_ro.setCaloShowerStep(remap.at(caloShowerStepPtr));

  _newCaloShowerStepROs->push_back(new_calo_shower_step_ro);
}

art::Ptr<mu2e::StepPointMC> mu2e::CompressDigiMCs::copyStepPointMC(const mu2e::StepPointMC& old_step, const InstanceLabel& instance) {

  keepSimParticle(old_step.simParticle());
  
  StepPointMC new_step(old_step);
  _newStepPointMCs.at(instance)->push_back(new_step);
  
  return art::Ptr<StepPointMC>(_newStepPointMCsPID.at(instance), _newStepPointMCs.at(instance)->size()-1, _newStepPointMCGetter.at(instance));
}

void mu2e::CompressDigiMCs::keepSimParticle(const art::Ptr<SimParticle>& sim_ptr) {

  // Also need to add all the parents too
  _simParticlesToKeep[sim_ptr.id()].insert(sim_ptr);
  art::Ptr<SimParticle> childPtr = sim_ptr;
  art::Ptr<SimParticle> parentPtr = childPtr->parent();
  
  while (parentPtr) {
    _simParticlesToKeep[sim_ptr.id()].insert(parentPtr);
    childPtr = parentPtr;
    parentPtr = parentPtr->parent();
  }
}


DEFINE_ART_MODULE(mu2e::CompressDigiMCs)
