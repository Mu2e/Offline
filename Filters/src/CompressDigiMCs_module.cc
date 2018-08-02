////////////////////////////////////////////////////////////////////////
// Class:       CompressDigiMCs
// Plugin Type: producer (art v2_06_02)
// File:        CompressDigiMCs_module.cc
//
// Creates new StrawDigiMC and CrvDigiMC collections after creating new
// StepPointMC, SimParticle, GenParticle and SimParticleTimeMaps with all 
// unnecessary MC objects removed
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

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "Mu2eUtilities/inc/compressSimParticleCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"

namespace mu2e {
  class CompressDigiMCs;

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
  art::Ptr<StepPointMC> copyStepPointMC(const mu2e::StepPointMC& old_step);

private:

  // art tags for the input collections
  art::InputTag _strawDigiMCTag;
  art::InputTag _crvDigiMCTag;

  std::vector<art::InputTag> _simParticleTags;
  std::vector<art::InputTag> _extraStepPointMCTags;
  std::vector<art::InputTag> _timeMapTags;
  art::InputTag _primarySimPtrsTag;

  // handles to the old collections
  art::Handle<StrawDigiMCCollection> _strawDigiMCsHandle;
  art::Handle<CrvDigiMCCollection> _crvDigiMCsHandle;
  std::vector<SimParticleTimeMap> _oldTimeMaps;
  art::Handle<SimParticlePtrCollection> _primarySimPtrsHandle;

  // unique_ptrs to the new output collections
  std::unique_ptr<StrawDigiMCCollection> _newStrawDigiMCs;
  std::unique_ptr<CrvDigiMCCollection> _newCrvDigiMCs;
  std::unique_ptr<StepPointMCCollection> _newStepPointMCs;
  std::map<art::ProductID, std::unique_ptr<SimParticleCollection> > _newSimParticles;
  std::map<art::ProductID, std::unique_ptr<GenParticleCollection> > _newGenParticles;
  std::vector<std::unique_ptr<SimParticleTimeMap> > _newSimParticleTimeMaps;
  std::unique_ptr<SimParticlePtrCollection> _newPrimarySimPtrs;

  // for StepPointMCs, SimParticles and GenParticles we also need reference their new locations with art::Ptrs and so need their ProductIDs and Getters
  art::ProductID _newStepPointMCsPID;
  const art::EDProductGetter* _newStepPointMCGetter;
  std::map<art::ProductID, art::ProductID> _newSimParticlesPID;
  std::map<art::ProductID, const art::EDProductGetter*> _oldSimParticleGetter;
  std::map<art::ProductID, const art::EDProductGetter*> _newSimParticleGetter;
  std::map<art::ProductID, art::ProductID> _newGenParticlesPID;
  std::map<art::ProductID, const art::EDProductGetter*> _newGenParticleGetter;

  // record the SimParticles that we are keeping so we can use compressSimParticleCollection to do all the work for us
  std::map<art::ProductID, SimParticleSelector> _simParticlesToKeep;
};


mu2e::CompressDigiMCs::CompressDigiMCs(fhicl::ParameterSet const & pset)
  : _strawDigiMCTag(pset.get<art::InputTag>("strawDigiMCTag")),
    _crvDigiMCTag(pset.get<art::InputTag>("crvDigiMCTag")),
    _simParticleTags(pset.get<std::vector<art::InputTag> >("simParticleTags")),
    _extraStepPointMCTags(pset.get<std::vector<art::InputTag> >("extraStepPointMCTags")),
    _timeMapTags(pset.get<std::vector<art::InputTag> >("timeMapTags")),
    _primarySimPtrsTag(pset.get<art::InputTag>("primarySimPtrsTag"))
{
  // Call appropriate produces<>() functions here.
  produces<StrawDigiMCCollection>();
  produces<CrvDigiMCCollection>();

  produces<StepPointMCCollection>();

  for (std::vector<art::InputTag>::const_iterator i_tag = _simParticleTags.begin(); i_tag != _simParticleTags.end(); ++i_tag) {
    produces<SimParticleCollection>( (*i_tag).label() );
    produces<GenParticleCollection>( (*i_tag).label() );
  }

  for (std::vector<art::InputTag>::const_iterator i_tag = _timeMapTags.begin(); i_tag != _timeMapTags.end(); ++i_tag) {
    produces<SimParticleTimeMap>( (*i_tag).label() );
  }

  produces<SimParticlePtrCollection>();
}

void mu2e::CompressDigiMCs::produce(art::Event & event)
{
  // Implementation of required member function here.
  _newStrawDigiMCs = std::unique_ptr<StrawDigiMCCollection>(new StrawDigiMCCollection);  
  _newCrvDigiMCs = std::unique_ptr<CrvDigiMCCollection>(new CrvDigiMCCollection);  
  _newStepPointMCs = std::unique_ptr<StepPointMCCollection>(new StepPointMCCollection);
  _newStepPointMCsPID = getProductID<StepPointMCCollection>();
  _newStepPointMCGetter = event.productGetter(_newStepPointMCsPID);

  // Create all the new collections, ProductIDs and product getters for the SimParticles and GenParticles
  // There is one for each background frame plus one for the primary event
  for (std::vector<art::InputTag>::const_iterator i_tag = _simParticleTags.begin(); i_tag != _simParticleTags.end(); ++i_tag) {
    const auto& oldSimParticles = event.getValidHandle<SimParticleCollection>(*i_tag);
    art::ProductID i_product_id = oldSimParticles.id();
    
    _simParticlesToKeep[i_product_id].clear();
    
    _newSimParticles[i_product_id] = std::unique_ptr<SimParticleCollection>(new SimParticleCollection);
    _newSimParticlesPID[i_product_id] = getProductID<SimParticleCollection>((*i_tag).label() );
    _newSimParticleGetter[i_product_id] = event.productGetter(_newSimParticlesPID[i_product_id]);

    _oldSimParticleGetter[i_product_id] = event.productGetter(i_product_id);
    
    _newGenParticles[i_product_id] = std::unique_ptr<GenParticleCollection>(new GenParticleCollection);
    _newGenParticlesPID[i_product_id] = getProductID<GenParticleCollection>((*i_tag).label() );
    _newGenParticleGetter[i_product_id] = event.productGetter(_newGenParticlesPID[i_product_id]);
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

  event.getByLabel(_primarySimPtrsTag, _primarySimPtrsHandle);
  const auto& primarySimPtrs = *_primarySimPtrsHandle;
  _newPrimarySimPtrs = std::unique_ptr<SimParticlePtrCollection>(new SimParticlePtrCollection);  


  event.getByLabel(_strawDigiMCTag, _strawDigiMCsHandle);
  const auto& strawDigiMCs = *_strawDigiMCsHandle;
  for (const auto& i_strawDigiMC : strawDigiMCs) {
    copyStrawDigiMC(i_strawDigiMC);
  }

  event.getByLabel(_crvDigiMCTag, _crvDigiMCsHandle);
  const auto& crvDigiMCs = *_crvDigiMCsHandle;
  for (const auto& i_crvDigiMC : crvDigiMCs) {
    copyCrvDigiMC(i_crvDigiMC);
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
	const SimParticleSelector& simParticles = simPartsToKeep.second;
	const std::set<cet::map_vector_key>& alreadyKeptKeys = simParticles.keys();
	for (const auto& alreadyKeptKey : alreadyKeptKeys) {
	  if (stepPointMC.simParticle()->id() == alreadyKeptKey) {
	    copyStepPointMC(stepPointMC);
	  }
	}
      }
    }
  }
  
  // Now compress the SimParticleCollections into their new collections
  for (std::vector<art::InputTag>::const_iterator i_tag = _simParticleTags.begin(); i_tag != _simParticleTags.end(); ++i_tag) {
    const auto& oldSimParticles = event.getValidHandle<SimParticleCollection>(*i_tag);
    art::ProductID i_product_id = oldSimParticles.id();
    compressSimParticleCollection(_newSimParticlesPID[i_product_id], _newSimParticleGetter[i_product_id], *oldSimParticles, 
				  _simParticlesToKeep[i_product_id], *(_newSimParticles[i_product_id]));

    for(auto& i : *(_newSimParticles[i_product_id])) {
      
      mu2e::SimParticle& newsim = i.second;
      if(!newsim.genParticle().isNull()) { // will crash if not resolvable

	// Copy GenParticle to the new collection
	_newGenParticles[i_product_id]->emplace_back(*newsim.genParticle());
	newsim.genParticle() = art::Ptr<GenParticle>(_newGenParticlesPID[i_product_id], _newGenParticles[i_product_id]->size()-1, _newGenParticleGetter[i_product_id]);
      }
      
      // Update the time maps
      art::Ptr<SimParticle> oldSimPtr(i_product_id, newsim.id().asUint(), _oldSimParticleGetter[i_product_id]);
      art::Ptr<SimParticle> newSimPtr(_newSimParticlesPID[i_product_id], newsim.id().asUint(), _newSimParticleGetter[i_product_id]);
      for (std::vector<SimParticleTimeMap>::const_iterator i_time_map = _oldTimeMaps.begin(); i_time_map != _oldTimeMaps.end(); ++i_time_map) {
	size_t i_element = i_time_map - _oldTimeMaps.begin();
	
	const SimParticleTimeMap& i_oldTimeMap = *i_time_map;
	SimParticleTimeMap& i_newTimeMap = *_newSimParticleTimeMaps.at(i_element);
	auto it = i_oldTimeMap.find(oldSimPtr);
	if (it != i_oldTimeMap.end()) {
	  i_newTimeMap[newSimPtr] = it->second;
	}
      }

      // Update the PrimarySimPtrs
      for (const auto& i_primarySimPtr : primarySimPtrs) {
	if (i_primarySimPtr == oldSimPtr) {
	  _newPrimarySimPtrs->push_back(newSimPtr);
	}
      }
    }
  }

  // Now add everything to the event
  event.put(std::move(_newStepPointMCs));  
  event.put(std::move(_newStrawDigiMCs));
  event.put(std::move(_newCrvDigiMCs));

  for (std::vector<art::InputTag>::const_iterator i_tag = _simParticleTags.begin(); i_tag != _simParticleTags.end(); ++i_tag) {
    const auto& oldSimParticles = event.getValidHandle<SimParticleCollection>(*i_tag);
    art::ProductID i_product_id = oldSimParticles.id();
    event.put(std::move(_newSimParticles[i_product_id]), (*i_tag).label());
    event.put(std::move(_newGenParticles[i_product_id]), (*i_tag).label());
  }

  for (std::vector<art::InputTag>::const_iterator i_tag = _timeMapTags.begin(); i_tag != _timeMapTags.end(); ++i_tag) {
    size_t i_element = i_tag - _timeMapTags.begin();
    event.put(std::move(_newSimParticleTimeMaps.at(i_element)), (*i_tag).label());
  }
  event.put(std::move(_newPrimarySimPtrs));
}

void mu2e::CompressDigiMCs::copyStrawDigiMC(const mu2e::StrawDigiMC& old_straw_digi_mc) {

  // Need to update the Ptrs for the StepPointMCs
  art::Ptr<StepPointMC> newTriggerStepPtr[StrawEnd::nends];
  for(int i_end=0;i_end<StrawEnd::nends;++i_end){
    StrawEnd::End end = static_cast<StrawEnd::End>(i_end);

    const art::Ptr<StepPointMC>& old_step_point = old_straw_digi_mc.stepPointMC(end);
    if (old_step_point.isAvailable()) {
      newTriggerStepPtr[i_end] = copyStepPointMC( *old_step_point );
    }
  }
  
  std::vector<art::Ptr<StepPointMC> > newWaveformStepPtrs;
  for (const auto& i_step_mc : old_straw_digi_mc.stepPointMCs()) {
    if (i_step_mc.isAvailable()) {
      newWaveformStepPtrs.push_back(copyStepPointMC(*i_step_mc));
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
      newStepPtrs.push_back(copyStepPointMC(*i_step_mc));
    }
  }

  art::Ptr<SimParticle> oldSimPtr = old_crv_digi_mc.GetSimParticle();
  art::Ptr<SimParticle> newSimPtr;
  if (oldSimPtr.isNonnull()) { // if the old CrvDigiMC doesn't have a null ptr for the SimParticle...
    newSimPtr = art::Ptr<SimParticle>(_newSimParticlesPID[oldSimPtr.id()], oldSimPtr->id().asUint(), _newSimParticleGetter[oldSimPtr.id()]);
  }
  else {
    newSimPtr = art::Ptr<SimParticle>();
  }
  CrvDigiMC new_crv_digi_mc(old_crv_digi_mc.GetVoltages(), newStepPtrs, 
			    newSimPtr, old_crv_digi_mc.GetStartTime(), 
			    old_crv_digi_mc.GetScintillatorBarIndex(), old_crv_digi_mc.GetSiPMNumber());

  _newCrvDigiMCs->push_back(new_crv_digi_mc);
}

art::Ptr<mu2e::StepPointMC> mu2e::CompressDigiMCs::copyStepPointMC(const mu2e::StepPointMC& old_step) {

  _simParticlesToKeep[old_step.simParticle().id()].push_back(old_step.simParticle()->id());
  art::Ptr<SimParticle> newSimPtr(_newSimParticlesPID[old_step.simParticle().id()], old_step.simParticle()->id().asUint(), _newSimParticleGetter[old_step.simParticle().id()]);

  // Also need to add all the parents (and change their genParticles) too
  art::Ptr<SimParticle> childPtr = old_step.simParticle();
  art::Ptr<SimParticle> parentPtr = childPtr->parent();
  
  while (parentPtr) {
    _simParticlesToKeep[old_step.simParticle().id()].push_back(parentPtr->id());

    childPtr = parentPtr;
    parentPtr = parentPtr->parent();
  }
  
  StepPointMC new_step(old_step);
  new_step.simParticle() = newSimPtr;

  _newStepPointMCs->push_back(new_step);
  
  return art::Ptr<StepPointMC>(_newStepPointMCsPID, _newStepPointMCs->size()-1, _newStepPointMCGetter);
}



DEFINE_ART_MODULE(mu2e::CompressDigiMCs)
