////////////////////////////////////////////////////////////////////////
// Class:       FilterKalSeed
// Plugin Type: filter (art v2_06_02)
// File:        FilterKalSeed_module.cc
//
// Generated at Wed Apr 12 16:10:46 2017 by Andrew Edmonds using cetskelgen
// from cetlib version v2_02_00.
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

#include <memory>

#include "DataProducts/inc/IndexMap.hh"
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/TrkQual.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "Mu2eUtilities/inc/compressSimParticleCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

namespace mu2e {
  class FilterKalSeed;

  class SimParticleSelector {
  public:
    SimParticleSelector() { }
    
    void push_back(cet::map_vector_key key) {
      m_keys.insert(key);
    }

    bool operator[]( cet::map_vector_key key ) const {
      return m_keys.find(key) != m_keys.end();
    }

  private:
    std::set<cet::map_vector_key> m_keys;
    
  };
}


class mu2e::FilterKalSeed : public art::EDFilter {
public:
  explicit FilterKalSeed(fhicl::ParameterSet const & pset);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FilterKalSeed(FilterKalSeed const &) = delete;
  FilterKalSeed(FilterKalSeed &&) = delete;
  FilterKalSeed & operator = (FilterKalSeed const &) = delete;
  FilterKalSeed & operator = (FilterKalSeed &&) = delete;

  // Required functions.
  bool filter(art::Event & event) override;

private:

  // Declare member data here.
  art::InputTag _kalFinalsTag;

  std::vector<art::InputTag> _simParticleTags;

  // unique_pts to the output collections
  std::unique_ptr<KalSeedCollection> _newKalFinals;
  std::unique_ptr<StrawDigiMCCollection> _newStrawDigiMCs;
  std::unique_ptr<StrawHitCollection> _newStrawHits;
  std::unique_ptr<StrawDigiCollection> _newStrawDigis;

  // for StepPointMCs and SimParticles we also need to copy and move objects to new collections
  std::unique_ptr<StepPointMCCollection> _newStepPointMCs;
  art::ProductID _newStepPointMCsPID;
  const art::EDProductGetter* _newStepPointMCGetter;

  // Need a different SimParticleCollection for each input collection to avoid two SimParticles with the same G4 TrackID being overwritten
  // These are mapped to the art::ProductID of the old SimParticleCollection
  std::map<art::ProductID, std::unique_ptr<SimParticleCollection> > _newSimParticles;
  std::map<art::ProductID, art::ProductID> _newSimParticlesPID;
  std::map<art::ProductID, const art::EDProductGetter*> _newSimParticleGetter;
  std::map<art::ProductID, SimParticleSelector> _simParticlesToKeep;

  // Need a different GenParticleCollection for each input collection to avoid two GenParticles with the same G4 TrackID being overwritten
  // These are mapped to the art::ProductID of the old SimParticleCollection (not GenParticleCollection)
  std::map<art::ProductID, std::unique_ptr<GenParticleCollection> > _newGenParticles;
  std::map<art::ProductID, art::ProductID> _newGenParticlesPID;
  std::map<art::ProductID, const art::EDProductGetter*> _newGenParticleGetter;

  void copyStrawDigiMC(const mu2e::StrawDigiMC& old_straw_digi_mc) {

    // Get information from the old StrawDigiMC
    double wetime[2] = {old_straw_digi_mc.wireEndTime(StrawDigi::TDCChannel::zero), old_straw_digi_mc.wireEndTime(StrawDigi::TDCChannel::one)};
    CLHEP::HepLorentzVector cpos[2] = {old_straw_digi_mc.clusterPosition(StrawDigi::TDCChannel::zero), old_straw_digi_mc.clusterPosition(StrawDigi::TDCChannel::one)};
    art::Ptr<StepPointMC> stepMC[2] = {old_straw_digi_mc.stepPointMC(StrawDigi::TDCChannel::zero), old_straw_digi_mc.stepPointMC(StrawDigi::TDCChannel::one)};

    art::Ptr<StepPointMC> newStepMCPtrs[2];
    for (int i_tdc = 0; i_tdc<2; ++i_tdc) {
      if(stepMC[i_tdc].isAvailable()) {
	newStepMCPtrs[i_tdc] = copyStepPointMC(*stepMC[i_tdc]);
      }
    }
    std::vector<art::Ptr<StepPointMC> > newStepMCs;
    for (const auto& i_step_mc : old_straw_digi_mc.stepPointMCs()) {
      newStepMCs.push_back(copyStepPointMC(*i_step_mc));
    }

    StrawDigiMC new_straw_digi_mc(old_straw_digi_mc.strawIndex(),
				  wetime,
				  cpos,
				  newStepMCPtrs,
				  newStepMCs);

    _newStrawDigiMCs->push_back(new_straw_digi_mc);
  }

  // Need to copy the StepPointMC to the new collection and update the art::Ptrs
  art::Ptr<StepPointMC> copyStepPointMC(const mu2e::StepPointMC& old_step) {

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

    StepPointMC new_step(newSimPtr,
			 old_step.volumeId(),
			 old_step.totalEDep(),
			 old_step.nonIonizingEDep(),
			 old_step.time(),
			 old_step.properTime(),
			 old_step.position(),
			 old_step.momentum(),
			 old_step.stepLength(),
			 old_step.endProcessCode());

    //    std::cout << old_step.simParticle()->id() << ": creationCode = " << old_step.simParticle()->creationCode() << ", PdgId = " << old_step.simParticle()->pdgId() << std::endl;
    _newStepPointMCs->push_back(new_step);

    return art::Ptr<StepPointMC>(_newStepPointMCsPID, _newStepPointMCs->size()-1, _newStepPointMCGetter);
  }

};


mu2e::FilterKalSeed::FilterKalSeed(fhicl::ParameterSet const & pset)
  : _kalFinalsTag(pset.get<art::InputTag>("kalFinalsTag", "KFFDeM"))
{
  // Call appropriate produces<>() functions here.
  produces<KalSeedCollection>();
  produces<StrawHitCollection>();
  produces<StrawDigiCollection>();
  produces<StrawDigiMCCollection>();
  produces<StepPointMCCollection>();
  produces<TrkQualCollection>();
  produces<IndexMap>();

  _simParticleTags.push_back("detectorFilter");
  _simParticleTags.push_back("protonMixer");
  _simParticleTags.push_back("deuteronMixer");
  _simParticleTags.push_back("photonMixer");
  _simParticleTags.push_back("neutronMixer");
  _simParticleTags.push_back("ootMixer");
  _simParticleTags.push_back("dioMixer");
  _simParticleTags.push_back("flashMixer");

  for (std::vector<art::InputTag>::const_iterator i_tag = _simParticleTags.begin(); i_tag != _simParticleTags.end(); ++i_tag) {
    produces<SimParticleCollection>( (*i_tag).label() );
    produces<GenParticleCollection>( (*i_tag).label() );
  }
}

bool mu2e::FilterKalSeed::filter(art::Event & event)
{
  // Implementation of required member function here.
  const auto& kalFinals = event.getValidHandle<KalSeedCollection>(_kalFinalsTag);
  const auto& strawHits = event.getValidHandle<StrawHitCollection>("makeSH");
  const auto& strawDigis = event.getValidHandle<StrawDigiCollection>("makeSD");
  const auto& strawDigiMCs = event.getValidHandle<StrawDigiMCCollection>("makeSD");
  
  _newKalFinals = std::unique_ptr<KalSeedCollection>(new KalSeedCollection);
  _newStrawHits = std::unique_ptr<StrawHitCollection>(new StrawHitCollection);
  _newStrawDigis = std::unique_ptr<StrawDigiCollection>(new StrawDigiCollection);
  _newStrawDigiMCs = std::unique_ptr<StrawDigiMCCollection>(new StrawDigiMCCollection);

  _newStepPointMCs = std::unique_ptr<StepPointMCCollection>(new StepPointMCCollection);
  _newStepPointMCsPID = getProductID<StepPointMCCollection>(event);
  _newStepPointMCGetter = event.productGetter(_newStepPointMCsPID);

  // Create all the new collections, ProductIDs and product getters for the SimParticles and GenParticles
  for (std::vector<art::InputTag>::const_iterator i_tag = _simParticleTags.begin(); i_tag != _simParticleTags.end(); ++i_tag) {
    const auto& oldSimParticles = event.getValidHandle<SimParticleCollection>(*i_tag);
    art::ProductID i_product_id = oldSimParticles.id();

    _newSimParticles[i_product_id] = std::unique_ptr<SimParticleCollection>(new SimParticleCollection);
    _newSimParticlesPID[i_product_id] = getProductID<SimParticleCollection>(event, (*i_tag).label() );
    _newSimParticleGetter[i_product_id] = event.productGetter(_newSimParticlesPID[i_product_id]);

    _newGenParticles[i_product_id] = std::unique_ptr<GenParticleCollection>(new GenParticleCollection);
    _newGenParticlesPID[i_product_id] = getProductID<GenParticleCollection>(event, (*i_tag).label() );
    _newGenParticleGetter[i_product_id] = event.productGetter(_newGenParticlesPID[i_product_id]);
  }

  std::unique_ptr<IndexMap> indexMap(new IndexMap);

  // Loop through the final kalman fits
  for (const auto& kalFinal : *kalFinals) {

    mu2e::CondensedIndex counter = 0;
    // Loop through the hits in the fit
    for (const auto& hit : kalFinal.hits()) {
      
      // at the moment, we are only writing everything related to the hits in the fit
      // in the future, we could add code here to look for other hits near these and write those out too
      StrawHitIndex hit_index = hit.index();

      mu2e::StrawHit straw_hit = (*strawHits).at(hit_index);
      _newStrawHits->push_back(straw_hit);

      mu2e::StrawDigi straw_digi = (*strawDigis).at(hit_index);
      _newStrawDigis->push_back(straw_digi);

      // Need a deep copy of StrawDigiMC
      const mu2e::StrawDigiMC& old_straw_digi_mc = (*strawDigiMCs).at(hit_index);
      copyStrawDigiMC(old_straw_digi_mc);

      // Update the IndexMap
      indexMap->addElement(hit.index(), counter);
      ++counter;
    }
    _newKalFinals->push_back(kalFinal); // want to keep all kal seeds
  }
  //  std::cout << *indexMap << std::endl;

  // Now compress the SimParticleCollections into their new collections
  for (std::vector<art::InputTag>::const_iterator i_tag = _simParticleTags.begin(); i_tag != _simParticleTags.end(); ++i_tag) {
    const auto& oldSimParticles = event.getValidHandle<SimParticleCollection>(*i_tag);
    art::ProductID i_product_id = oldSimParticles.id();

    compressSimParticleCollection(_newSimParticlesPID[i_product_id], _newSimParticleGetter[i_product_id], *oldSimParticles, _simParticlesToKeep[i_product_id], *(_newSimParticles[i_product_id]));

    for(auto& i : *(_newSimParticles[i_product_id])) {
      mu2e::SimParticle& newsim = i.second;
      if(!newsim.genParticle().isNull()) { // will crash if not resolvable

	// Copy GenParticle to the new collection
	_newGenParticles[i_product_id]->emplace_back(*newsim.genParticle());
	newsim.genParticle() = art::Ptr<GenParticle>(_newGenParticlesPID[i_product_id], _newGenParticles[i_product_id]->size()-1, _newGenParticleGetter[i_product_id]);
      }
    }
  }
  event.put(std::move(_newKalFinals));
  event.put(std::move(_newStrawHits));
  event.put(std::move(_newStrawDigis));
  event.put(std::move(_newStrawDigiMCs));

  for (std::vector<art::InputTag>::const_iterator i_tag = _simParticleTags.begin(); i_tag != _simParticleTags.end(); ++i_tag) {
    const auto& oldSimParticles = event.getValidHandle<SimParticleCollection>(*i_tag);
    art::ProductID i_product_id = oldSimParticles.id();
    event.put(std::move(_newSimParticles[i_product_id]), (*i_tag).label());
    event.put(std::move(_newGenParticles[i_product_id]), (*i_tag).label());
  }
  event.put(std::move(indexMap));


  // Get the hits from the virtualdetector
  const auto& stepPointMCs = event.getValidHandle<StepPointMCCollection>("detectorFilter:virtualdetector");
  for (const auto& stepPointMC : *stepPointMCs) {
    if (stepPointMC.volumeId() == VirtualDetectorId::TT_FrontHollow) { // can add other virtualdetectors here
      copyStepPointMC(stepPointMC);
    }
  }
  event.put(std::move(_newStepPointMCs));
  

  // Get the TrkQual information that was used for these KalSeeds
  const auto& trkQuals = event.getValidHandle<TrkQualCollection>("KFFDeM");
  std::unique_ptr<TrkQualCollection> newTrkQuals(new TrkQualCollection);
  for (const auto& trkQual : *trkQuals) {
    newTrkQuals->push_back(trkQual);
  }

  // Put everything into the event
  event.put(std::move(newTrkQuals));
  

  return true;
}

DEFINE_ART_MODULE(mu2e::FilterKalSeed)
