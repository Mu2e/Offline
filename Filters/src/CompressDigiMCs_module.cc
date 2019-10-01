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
#include "art_root_io/TFileService.h"

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
#include "DataProducts/inc/IndexMap.hh"
#include "MCDataProducts/inc/CaloClusterMC.hh"
#include "MCDataProducts/inc/CrvCoincidenceClusterMCCollection.hh"
#include "MCDataProducts/inc/PrimaryParticle.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"

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
  typedef std::map<art::Ptr<mu2e::StepPointMC>, art::Ptr<mu2e::StepPointMC> > StepPointMCRemap;
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
  void copyCaloClusterMC(const mu2e::CaloClusterMC& old_calo_cluster_mc);
  void copyCrvCoincClusterMC(const mu2e::CrvCoincidenceClusterMC& old_crv_coinc_cluster_mc);
  void copyPrimaryParticle(const mu2e::PrimaryParticle& old_primary_particle);

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

  // Optional parameters for reco output
  art::InputTag _strawDigiMCIndexMapTag;
  mu2e::IndexMap _strawDigiMCIndexMap;
  art::InputTag _crvDigiMCIndexMapTag;
  mu2e::IndexMap _crvDigiMCIndexMap;
  art::InputTag _caloClusterMCTag;
  art::Handle<CaloClusterMCCollection> _caloClusterMCsHandle;
  std::unique_ptr<CaloClusterMCCollection> _newCaloClusterMCs;
  bool _keepAllGenParticles;
  art::InputTag _crvCoincClusterMCTag;
  art::Handle<CrvCoincidenceClusterMCCollection> _crvCoincClusterMCsHandle;
  std::unique_ptr<CrvCoincidenceClusterMCCollection> _newCrvCoincClusterMCs;
  art::InputTag _primaryParticleTag;
  art::Handle<PrimaryParticle> _primaryParticleHandle;
  std::unique_ptr<PrimaryParticle> _newPrimaryParticle;
  bool _rekeySimParticleCollection;

  // other optional parameters
  art::InputTag _mcTrajectoryTag;
  art::Handle<MCTrajectoryCollection> _mcTrajectoriesHandle;
  std::unique_ptr<MCTrajectoryCollection> _newMCTrajectories;

  // For CrvDigiMCs, there's a chance that the same StepPointMC will go into multiple CrvDigiMCs
  // This module didn't take this into account initially and so the same StepPointMC was being written out multiple times
  // This std::set is used to make sure that this doesn't happen
  std::set<art::Ptr<StepPointMC> > _crvStepPointMCsSeen;
  std::map<art::Ptr<StepPointMC>, art::Ptr<StepPointMC> > _crvStepPointMCsMap;
};


mu2e::CompressDigiMCs::CompressDigiMCs(fhicl::ParameterSet const & pset)
  : art::EDProducer{pset},
    _strawDigiMCTag(pset.get<art::InputTag>("strawDigiMCTag")),
    _crvDigiMCTag(pset.get<art::InputTag>("crvDigiMCTag")),
    _simParticleTags(pset.get<std::vector<art::InputTag> >("simParticleTags")),
    _extraStepPointMCTags(pset.get<std::vector<art::InputTag> >("extraStepPointMCTags")),
    _timeMapTags(pset.get<std::vector<art::InputTag> >("timeMapTags")),
    _caloShowerStepTags(pset.get<std::vector<art::InputTag> >("caloShowerStepTags")),
    _caloShowerSimTag(pset.get<art::InputTag>("caloShowerSimTag")),
    _caloShowerStepROTag(pset.get<art::InputTag>("caloShowerStepROTag")),
    _trackerOutputInstanceLabel(pset.get<std::string>("trackerOutputInstanceLabel", "tracker")),
    _crvOutputInstanceLabel(pset.get<std::string>("crvOutputInstanceLabel", "CRV")),
    _strawDigiMCIndexMapTag(pset.get<art::InputTag>("strawDigiMCIndexMapTag", "")),
    _crvDigiMCIndexMapTag(pset.get<art::InputTag>("crvDigiMCIndexMapTag", "")),
    _caloClusterMCTag(pset.get<art::InputTag>("caloClusterMCTag", "")),
    _keepAllGenParticles(pset.get<bool>("keepAllGenParticles", true)),
    _crvCoincClusterMCTag(pset.get<art::InputTag>("crvCoincClusterMCTag", "")),
    _primaryParticleTag(pset.get<art::InputTag>("primaryParticleTag", "")),
    _rekeySimParticleCollection(pset.get<bool>("rekeySimParticleCollection", true)),
    _mcTrajectoryTag(pset.get<art::InputTag>("mcTrajectoryTag", ""))
{
  // Call appropriate produces<>() functions here.
  produces<StrawDigiMCCollection>();
  produces<CrvDigiMCCollection>();

  _newStepPointMCInstances.push_back(_trackerOutputInstanceLabel); // will always be a tracker and CRV instance
  _newStepPointMCInstances.push_back(_crvOutputInstanceLabel);     // filled with the StepPointMCs referenced by their DigiMCs
  for (std::vector<art::InputTag>::const_iterator i_tag = _extraStepPointMCTags.begin(); i_tag != _extraStepPointMCTags.end(); ++i_tag) {
    _newStepPointMCInstances.push_back( (*i_tag).instance() );
  }

  for (const auto& i_instance : _newStepPointMCInstances) {
    produces<StepPointMCCollection>( i_instance );
  }

  produces<SimParticleCollection>();
  produces<GenParticleCollection>();

  for (std::vector<art::InputTag>::const_iterator i_tag = _timeMapTags.begin(); i_tag != _timeMapTags.end(); ++i_tag) {
    if ( (*i_tag).instance() == "") {
      produces<SimParticleTimeMap>( (*i_tag).label() ); // in the case where we have two time maps from different modules
    }
    else {
      produces<SimParticleTimeMap>( (*i_tag).instance() ); // in the case where we have two time maps from the same module (e.g. a first stage compression)
    }
  }

  // Two possible compressions for calorimeter
  if (_caloClusterMCTag == "") {
    produces<CaloShowerStepCollection>();
    produces<CaloShowerSimCollection>();
    produces<CaloShowerStepROCollection>();
  }
  if (_caloClusterMCTag != "") {
    produces<CaloClusterMCCollection>();
  }
  if (_crvCoincClusterMCTag != "") {
    produces<CrvCoincidenceClusterMCCollection>();
  }
  if (_primaryParticleTag != "") {
    produces<PrimaryParticle>();
  }
  if (_mcTrajectoryTag != "") {
    produces<MCTrajectoryCollection>();
  }
}

void mu2e::CompressDigiMCs::produce(art::Event & event)
{
  // Implementation of required member function here.
  _newStrawDigiMCs = std::unique_ptr<StrawDigiMCCollection>(new StrawDigiMCCollection);
  _newCrvDigiMCs = std::unique_ptr<CrvDigiMCCollection>(new CrvDigiMCCollection);

  for (const auto& i_instance : _newStepPointMCInstances) {
    _newStepPointMCs[i_instance] = std::unique_ptr<StepPointMCCollection>(new StepPointMCCollection);
    _newStepPointMCsPID[i_instance] = event.getProductID<StepPointMCCollection>(i_instance);
    _newStepPointMCGetter[i_instance] = event.productGetter(_newStepPointMCsPID[i_instance]);
  }

  _newSimParticles = std::unique_ptr<SimParticleCollection>(new SimParticleCollection);
  _newSimParticlesPID = event.getProductID<SimParticleCollection>();
  _newSimParticleGetter = event.productGetter(_newSimParticlesPID);

  _newGenParticles = std::unique_ptr<GenParticleCollection>(new GenParticleCollection);
  _newGenParticlesPID = event.getProductID<GenParticleCollection>();
  _newGenParticleGetter = event.productGetter(_newGenParticlesPID);

  // Create all the new collections, ProductIDs and product getters for the SimParticles and GenParticles
  // There is one for each background frame plus one for the primary event
  unsigned int n_gen_particles_to_keep = 0;
  for (std::vector<art::InputTag>::const_iterator i_tag = _simParticleTags.begin(); i_tag != _simParticleTags.end(); ++i_tag) {
    const auto& oldSimParticles = event.getValidHandle<SimParticleCollection>(*i_tag);
    art::ProductID i_product_id = oldSimParticles.id();
    const art::EDProductGetter* i_product_getter = event.productGetter(i_product_id);

    _simParticlesToKeep[i_product_id].clear();

    if (_keepAllGenParticles) {
      // Add all the SimParticles that are also GenParticles
      for (const auto& i_oldSimParticle : *oldSimParticles) {
        const cet::map_vector_key& key = i_oldSimParticle.first;
        const SimParticle& i_oldSim = i_oldSimParticle.second;
        if (i_oldSim.genParticle().isNonnull()) {
          keepSimParticle(art::Ptr<SimParticle>(i_product_id, key.asUint(), i_product_getter));
          ++n_gen_particles_to_keep;
        }
      }
    }
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

  // If we've been given IndexMap tags, then use those
  if (_strawDigiMCIndexMapTag != "") {
    art::Handle<mu2e::IndexMap> indexMapHandle;
    event.getByLabel(_strawDigiMCIndexMapTag, indexMapHandle);
    _strawDigiMCIndexMap = *indexMapHandle;
  }
  if (_crvDigiMCIndexMapTag != "") {
    art::Handle<mu2e::IndexMap> indexMapHandle;
    event.getByLabel(_crvDigiMCIndexMapTag, indexMapHandle);
    _crvDigiMCIndexMap = *indexMapHandle;
  }
  // If we have a CaloClusterMC collection, use that
  if (_caloClusterMCTag != "") {
    event.getByLabel(_caloClusterMCTag, _caloClusterMCsHandle);
    _newCaloClusterMCs = std::unique_ptr<CaloClusterMCCollection>(new CaloClusterMCCollection);
  }
  // If we have a CrvCoincClusterMC collection, use that
  if (_crvCoincClusterMCTag != "") {
    event.getByLabel(_crvCoincClusterMCTag, _crvCoincClusterMCsHandle);
    _newCrvCoincClusterMCs = std::unique_ptr<CrvCoincidenceClusterMCCollection>(new CrvCoincidenceClusterMCCollection);
  }
  // If we have a PrimaryParticle, use that
  if (_primaryParticleTag != "") {
    event.getByLabel(_primaryParticleTag, _primaryParticleHandle);
    _newPrimaryParticle = std::unique_ptr<PrimaryParticle>(new PrimaryParticle);
  }
  // If we want to keep MC trajectories
  if (_mcTrajectoryTag != "") {
    event.getByLabel(_mcTrajectoryTag, _mcTrajectoriesHandle);
    _newMCTrajectories = std::unique_ptr<MCTrajectoryCollection>(new MCTrajectoryCollection);
  }


  // Now start to compress
  event.getByLabel(_strawDigiMCTag, _strawDigiMCsHandle);
  const auto& strawDigiMCs = *_strawDigiMCsHandle;
  for (size_t i = 0; i < strawDigiMCs.size(); ++i) {
    const auto& i_strawDigiMC = strawDigiMCs.at(i);
    mu2e::FullIndex full_i = i;
    bool in_index_map = false;
    if (_strawDigiMCIndexMapTag != "") {
      in_index_map = _strawDigiMCIndexMap.checkInMap(full_i);
    }
    if (_strawDigiMCIndexMapTag == "" || in_index_map) {
      copyStrawDigiMC(i_strawDigiMC);
    }
  }

  // Only check for this if we are not reducing the number of StrawDigiMCs
  if (_strawDigiMCIndexMapTag == "" && strawDigiMCs.size() != _newStrawDigiMCs->size()) {
    throw cet::exception("CompressDigiMCs") << "The number of StrawDigiMCs before and after compression does not match ("
                                            << strawDigiMCs.size() << " != " << _newStrawDigiMCs->size() << ")" << std::endl;
  }


  if (_crvDigiMCTag != "") {
    _crvStepPointMCsSeen.clear();
    _crvStepPointMCsMap.clear();

    event.getByLabel(_crvDigiMCTag, _crvDigiMCsHandle);
    const auto& crvDigiMCs = *_crvDigiMCsHandle;
    for (size_t i = 0; i < crvDigiMCs.size(); ++i) {
      const auto& i_crvDigiMC = crvDigiMCs.at(i);
      mu2e::FullIndex full_i = i;
      bool in_index_map = false;
      if (_crvDigiMCIndexMapTag != "") {
        in_index_map = _crvDigiMCIndexMap.checkInMap(full_i);
      }
      if (_crvDigiMCIndexMapTag == "" || in_index_map) {
        copyCrvDigiMC(i_crvDigiMC);
      }
    }
    // Only check for this if we are not reducing the number of CrvDigiMCs
    if (_crvDigiMCIndexMapTag == "" && crvDigiMCs.size() != _newCrvDigiMCs->size()) {
      throw cet::exception("CompressDigiMCs") << "The number of CrvDigiMCs before and after compression does not match ("
                                              << crvDigiMCs.size() << " != " << _newCrvDigiMCs->size() << ")" << std::endl;
    }
  }

  // Two possible compressions for calorimeter
  // The first just takes the CaloShowerSteps, CaloShowerSims and CaloShowerStepROs and reassigns Ptrs (i.e. no actual compression....)
  if (_caloClusterMCTag == "") {
    CaloShowerStepRemap caloShowerStepRemap;
    _newCaloShowerSteps = std::unique_ptr<CaloShowerStepCollection>(new CaloShowerStepCollection);
    _newCaloShowerStepsPID = event.getProductID<CaloShowerStepCollection>();
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
  }
  // The second uses CaloClusterMCs and only keeps SimParticles that have been assigned
  else if (_caloClusterMCTag != "") {
    const auto& caloClusterMCs = *_caloClusterMCsHandle;
    for (const auto& i_caloClusterMC : caloClusterMCs) {
      copyCaloClusterMC(i_caloClusterMC);
    }
  }

  // Optional CrvCoincidenceClusterMCs
  if (_crvCoincClusterMCTag != "") {
    const auto& crvCoincClusterMCs = *_crvCoincClusterMCsHandle;
    for (const auto& i_crvCoincClusterMC : crvCoincClusterMCs) {
      copyCrvCoincClusterMC(i_crvCoincClusterMC);
    }
  }

  // Optional primary particles
  if (_primaryParticleTag != "") {
    const auto& primaryParticle = *_primaryParticleHandle;
    copyPrimaryParticle(primaryParticle);
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
    if (_rekeySimParticleCollection) {
      compressSimParticleCollection(_newSimParticlesPID, _newSimParticleGetter, *oldSimParticles,
                                    simPartSelector, *_newSimParticles, keyRemap);
    }
    else {
      compressSimParticleCollection(_newSimParticlesPID, _newSimParticleGetter, *oldSimParticles,
                                    simPartSelector, *_newSimParticles);
    }

    // Fill out the SimParticleRemapping
    for (const auto& i_keptSimPart : _simParticlesToKeep[i_product_id]) {
      cet::map_vector_key oldKey = cet::map_vector_key(i_keptSimPart.key());
      cet::map_vector_key newKey = oldKey;
      if (_rekeySimParticleCollection) {
        newKey = keyRemap->at(oldKey);
      }
      remap[i_keptSimPart] = art::Ptr<SimParticle>(_newSimParticlesPID, newKey.asUint(), _newSimParticleGetter);
    }
  }
  if (keep_size != _newSimParticles->size()) {
    throw cet::exception("CompressDigiMCs") << "Number of SimParticles in output collection ("
                                            << _newSimParticles->size()
                                            << ") does not match the number of SimParticles we wanted to keep ("
                                            << keep_size << ")" << std::endl;
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
  if (_keepAllGenParticles && _newGenParticles->size() != n_gen_particles_to_keep) {
    throw cet::exception("CompressDigiMCs") << "Number of GenParticles in output collection does not match the number of GenParticles we wanted to keep (" << n_gen_particles_to_keep << " != " << _newGenParticles->size() << ")" << std::endl;
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

  if (_caloClusterMCTag == "") {
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
  }
  else if (_caloClusterMCTag != "") {
    for (auto& i_caloClusterMC : *_newCaloClusterMCs) {
      for (auto& i_caloMCEDep : i_caloClusterMC._edeps) {
        art::Ptr<SimParticle> newSimPtr = remap.at(i_caloMCEDep.sim());
        i_caloMCEDep._simp = newSimPtr;
      }
    }
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
  // Update CrvCoincClusterMCs if needs be
  if (_crvCoincClusterMCTag != "") {
    for (auto& i_crvCoincClusterMC : *_newCrvCoincClusterMCs) {
      for (auto& i_pulseInfo : i_crvCoincClusterMC.GetModifiablePulses()) {
        art::Ptr<SimParticle> oldSimPtr = i_pulseInfo._simParticle;
        art::Ptr<SimParticle> newSimPtr = remap.at(oldSimPtr);
        i_pulseInfo._simParticle = newSimPtr;
      }

      art::Ptr<SimParticle> oldSimPtr = i_crvCoincClusterMC.GetMostLikelySimParticle();
      art::Ptr<SimParticle> newSimPtr = remap.at(oldSimPtr);
      i_crvCoincClusterMC.SetMostLikelySimParticle(newSimPtr);
    }
  }
  // Update PrimaryParticle if needs be
  if (_primaryParticleTag != "") {
    for (auto& i_simPartPtr : _newPrimaryParticle->modifySimParticles()) {
      i_simPartPtr = remap.at(i_simPartPtr);
    }
  }
  // Create new MC Trajectory collection
  if (_mcTrajectoryTag != "") {
    for (const auto& i_mcTrajectory : *_mcTrajectoriesHandle) {
      art::Ptr<SimParticle> oldSimPtr = i_mcTrajectory.first;
      if (remap.find(oldSimPtr) != remap.end()) {
        _newMCTrajectories->insert(std::pair<art::Ptr<SimParticle>, mu2e::MCTrajectory>(remap.at(oldSimPtr), i_mcTrajectory.second));
      }
    }
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
    if ( (*i_tag).instance() == "") {
      event.put(std::move(_newSimParticleTimeMaps.at(i_element)), (*i_tag).label());
    }
    else {
      event.put(std::move(_newSimParticleTimeMaps.at(i_element)), (*i_tag).instance());
    }
  }

  if (_caloClusterMCTag == "") {
    event.put(std::move(_newCaloShowerSteps));
    event.put(std::move(_newCaloShowerSims));
    event.put(std::move(_newCaloShowerStepROs));
  }
  else if (_caloClusterMCTag != "") {
    event.put(std::move(_newCaloClusterMCs));
  }

  if (_crvCoincClusterMCTag != "") {
    event.put(std::move(_newCrvCoincClusterMCs));
  }
  if (_primaryParticleTag != "") {
    event.put(std::move(_newPrimaryParticle));
  }
  if (_mcTrajectoryTag != "") {
    event.put(std::move(_newMCTrajectories));
  }
}

void mu2e::CompressDigiMCs::copyStrawDigiMC(const mu2e::StrawDigiMC& old_straw_digi_mc) {

  StepPointMCRemap step_remap;

  // Need to update the Ptrs for the StepPointMCs
  art::Ptr<StepPointMC> newTriggerStepPtr[StrawEnd::nends];
  for(int i_end=0;i_end<StrawEnd::nends;++i_end){
    StrawEnd::End end = static_cast<StrawEnd::End>(i_end);

    const art::Ptr<StepPointMC>& old_step_point = old_straw_digi_mc.stepPointMC(end);
    const auto& newStepPtrIter = step_remap.find(old_step_point);
    if (newStepPtrIter == step_remap.end()) {
      if (old_step_point.isAvailable()) {
	step_remap[old_step_point] = copyStepPointMC( *old_step_point, _trackerOutputInstanceLabel );
      }
      else { // this is a null Ptr but it should be added anyway to keep consistency (not expected for StrawDigis)
	step_remap[old_step_point] = old_step_point;
      }
    }
    art::Ptr<StepPointMC> new_step_point = step_remap.at(old_step_point);
    newTriggerStepPtr[i_end] = new_step_point;
  }

  std::vector<art::Ptr<StepPointMC> > newWaveformStepPtrs;
  for (const auto& i_step_mc : old_straw_digi_mc.stepPointMCs()) {
    const auto& newStepPtrIter = step_remap.find(i_step_mc);
    if (newStepPtrIter == step_remap.end()) {
      if (i_step_mc.isAvailable()) {
	step_remap[i_step_mc] = copyStepPointMC( *i_step_mc, _trackerOutputInstanceLabel );
      }
      else { // this is a null Ptr but it should be added anyway to keep consistency (not expected for StrawDigis)
	step_remap[i_step_mc] = i_step_mc;
      }
    }
    art::Ptr<StepPointMC> new_step_point = step_remap.at(i_step_mc);
    newWaveformStepPtrs.push_back(new_step_point);
  }

  StrawDigiMC new_straw_digi_mc(old_straw_digi_mc, newTriggerStepPtr, newWaveformStepPtrs); // copy everything except the Ptrs from the old StrawDigiMC
  _newStrawDigiMCs->push_back(new_straw_digi_mc);
}

void mu2e::CompressDigiMCs::copyCrvDigiMC(const mu2e::CrvDigiMC& old_crv_digi_mc) {

  // Need to update the Ptrs for the StepPointMCs
  std::vector<art::Ptr<StepPointMC> > newStepPtrs;
  for (const auto& i_step_mc : old_crv_digi_mc.GetStepPoints()) {
    if (i_step_mc.isAvailable()) {
      if (_crvStepPointMCsSeen.insert(i_step_mc).second == true) { // if we have inserted this StepPointMCPtrs (i.e. it hasn't already been seen)
        art::Ptr<StepPointMC> newStepPtr = copyStepPointMC(*i_step_mc, _crvOutputInstanceLabel);
        newStepPtrs.push_back(newStepPtr);
        _crvStepPointMCsMap[i_step_mc] = newStepPtr;
      }
      else {
        newStepPtrs.push_back(_crvStepPointMCsMap.at(i_step_mc));
      }
    }
    else { // this is a null Ptr but it should be added anyway to keep consistency (expected for CrvDigis)
      newStepPtrs.push_back(i_step_mc);
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

    return art::Ptr<mu2e::CaloShowerStep>(_newCaloShowerStepsPID, _newCaloShowerSteps->size()-1, _newCaloShowerStepGetter);
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

void mu2e::CompressDigiMCs::copyCaloClusterMC(const mu2e::CaloClusterMC& old_calo_cluster_mc) {

  for (const auto& i_caloMCEDep : old_calo_cluster_mc.energyDeposits()) {
    keepSimParticle(i_caloMCEDep.sim());
  }
  CaloClusterMC new_calo_cluster_mc(old_calo_cluster_mc);
  _newCaloClusterMCs->push_back(new_calo_cluster_mc);
}

void mu2e::CompressDigiMCs::copyCrvCoincClusterMC(const mu2e::CrvCoincidenceClusterMC& old_crv_coinc_cluster_mc) {

  for (const auto& i_pulseInfo : old_crv_coinc_cluster_mc.GetPulses()) {
    keepSimParticle(i_pulseInfo._simParticle);
  }
  keepSimParticle(old_crv_coinc_cluster_mc.GetMostLikelySimParticle());

  CrvCoincidenceClusterMC new_crv_coinc_cluster_mc(old_crv_coinc_cluster_mc);
  _newCrvCoincClusterMCs->push_back(new_crv_coinc_cluster_mc);
}

void mu2e::CompressDigiMCs::copyPrimaryParticle(const mu2e::PrimaryParticle& old_primary_particle) {

  for (const auto& i_simPart : old_primary_particle.primarySimParticles()) {
    keepSimParticle(i_simPart);
  }

  *_newPrimaryParticle = old_primary_particle;
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
