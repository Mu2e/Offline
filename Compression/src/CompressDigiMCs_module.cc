////////////////////////////////////////////////////////////////////////
// Class:       CompressDigiMCs
// Plugin Type: producer (art v2_06_02)
// File:        CompressDigiMCs_module.cc
//
// Creates new StrawDigiMC and CrvDigiMC collections after creating new
// StepPointMC, SimParticle, and GenParticle with all
// unnecessary MC objects removed.
//
// Also creates new CaloShowerStep, CaloShowerRO and CaloShowerSim collections after
// remapping the art::Ptrs to the SimParticles
//
// Generated at Wed Apr 12 16:10:46 2017 by Andrew Edmonds using cetskelgen
// from cetlib version v2_02_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"

#include <memory>

#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/CrvDigiMC.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/CaloShowerSim.hh"
#include "Offline/MCDataProducts/inc/CaloShowerRO.hh"
#include "Offline/MCDataProducts/inc/CaloClusterMC.hh"

#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/Mu2eUtilities/inc/compressSimParticleCollection.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticleRemapping.hh"
#include "Offline/DataProducts/inc/IndexMap.hh"
#include "Offline/MCDataProducts/inc/CrvCoincidenceClusterMC.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/MCDataProducts/inc/SurfaceStep.hh"

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

  typedef std::string InstanceLabel;
  typedef std::map<cet::map_vector_key, cet::map_vector_key> KeyRemap;
  typedef std::map<art::Ptr<mu2e::StepPointMC>, art::Ptr<mu2e::StepPointMC> > StepPointMCRemap;
  typedef std::map<art::Ptr<mu2e::StrawGasStep>, art::Ptr<mu2e::StrawGasStep> > StrawGasStepRemap;
  typedef std::map<art::Ptr<mu2e::CaloShowerStep>, art::Ptr<mu2e::CaloShowerStep> > CaloShowerStepRemap;
  typedef std::map<art::Ptr<mu2e::CrvStep>, art::Ptr<mu2e::CrvStep> > CrvStepRemap;
  typedef std::map<art::Ptr<mu2e::SurfaceStep>, art::Ptr<mu2e::SurfaceStep> > SurfaceStepRemap;
}


class mu2e::CompressDigiMCs : public art::EDProducer {
public:
  struct Config {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;

    // fhicl parameters for all of our inputs
    fhicl::Atom<art::InputTag> strawDigiMCTag{Name("strawDigiMCTag"), Comment("InputTag for the StrawDigiMCCollection")};
    fhicl::Atom<art::InputTag> crvDigiMCTag{Name("crvDigiMCTag"), Comment("InputTag for the CrvDigiMCCollection")};
    fhicl::Sequence<art::InputTag> simParticleTags{Name("simParticleTags"), Comment("Sequence of InputTags to the SimParticleCollections")};
    fhicl::Sequence<art::InputTag> extraStepPointMCTags{Name("extraStepPointMCTags"), Comment("Sequence of InputTags for additional StepPointMCCollections that you want to keep the steps from")};
    fhicl::Sequence<art::InputTag> surfaceStepTags{Name("surfaceStepTags"), Comment("InputTags for SurfaceSteps we want to keep")};
    fhicl::Sequence<art::InputTag> caloShowerStepTags{Name("caloShowerStepTags"), Comment("Sequence of InputTags for CaloShowerSteps")};
    fhicl::Atom<art::InputTag> caloShowerSimTag{Name("caloShowerSimTag"), Comment("InputTag for the CaloShowerSim")};
    fhicl::Atom<art::InputTag> caloShowerROTag{Name("caloShowerROTag"), Comment("InputTag for the CaloShowerRO")};

    fhicl::Atom<art::InputTag> strawDigiMCIndexMapTag{Name("strawDigiMCIndexMapTag"), Comment("InputTag for an IndexMap that maps the StrawDigiMCs we want to keep to those in the uncompressed StrawDigiMCCollection (leave blank to keep all StrawDigiMCs")};
    fhicl::Atom<art::InputTag> crvDigiMCIndexMapTag{Name("crvDigiMCIndexMapTag"), Comment("InputTag for an IndexMap that maps the CrvDigiMCs we want to keep to those in the uncompressed CrvDigiMCCollection (leave blank to keep all CrvDigiMCs")};

    // Reco objects
    fhicl::Atom<art::InputTag> caloClusterMCTag{Name("caloClusterMCTag"), Comment("InputTag for CaloClusterMCCollection")};
    fhicl::Sequence<art::InputTag> crvCoincClusterMCTags{Name("crvCoincClusterMCTags"), Comment("InputTags for CrvCoincidenceClusterMCCollections")};

    fhicl::Atom<art::InputTag> primaryParticleTag{Name("primaryParticleTag"), Comment("InputTag for PrimarParticle")};
    fhicl::Atom<art::InputTag> mcTrajectoryTag{Name("mcTrajectoryTag"), Comment("InputTag for the MCTrajectoryCollection")};

    // fhicl parameters for output
    fhicl::Atom<bool> keepAllGenParticles{Name("keepAllGenParticles"), Comment("Set to true if you want to keep all GenParticles even if their descendents make no hits in the detector")};
    fhicl::Atom<bool> rekeySimParticleCollection{Name("rekeySimParticleCollection"), Comment("Set to true to change the keys in the SimParticleCollection (necessary for mixed events)")};

    fhicl::Atom<bool> noCompression{Name("noCompression"), Comment("Set to true to turn off compression"), false};

    // detector steps we may want to keep all of
    fhicl::Sequence<art::InputTag> crvStepsToKeep{Name("crvStepsToKeep"), Comment("InputTags for CrvSteps we want to keep")};
  };
  typedef art::EDProducer::Table<Config> Parameters;

  explicit CompressDigiMCs(const Parameters& conf);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Required functions.
  void produce(art::Event & event) override;

  // Other functions
  void copyStrawDigiMC(const mu2e::StrawDigiMC& old_straw_digi_mc);
  void copyCrvDigiMC(const mu2e::CrvDigiMC& old_crv_digi_mc);
  art::Ptr<StepPointMC> copyStepPointMC(const mu2e::StepPointMC& old_step, const InstanceLabel& instance);
  art::Ptr<StrawGasStep> copyStrawGasStep(const mu2e::StrawGasStep& old_step);
  art::Ptr<CrvStep> copyCrvStep(const mu2e::CrvStep& old_step);
  art::Ptr<SurfaceStep> copySurfaceStep(const mu2e::SurfaceStep& old_step);
  art::Ptr<mu2e::CaloShowerStep> copyCaloShowerStep(const mu2e::CaloShowerStep& old_calo_shower_step);
  void copyCaloShowerSim(const mu2e::CaloShowerSim& old_calo_shower_sim, const CaloShowerStepRemap& remap);
  void copyCaloShowerRO(const mu2e::CaloShowerRO& old_calo_shower_step_ro, const CaloShowerStepRemap& remap);
  void keepSimParticle(const art::Ptr<SimParticle>& sim_ptr);
  void copyCaloClusterMC(const mu2e::CaloClusterMC& old_calo_cluster_mc);
  art::Ptr<CaloHitMC> copyCaloHitMC(const mu2e::CaloHitMC& old_calo_hit_mc);
  void copyCrvCoincClusterMC(const mu2e::CrvCoincidenceClusterMC& old_crv_coinc_cluster_mc, size_t i_tag);
  void copyPrimaryParticle(const mu2e::PrimaryParticle& old_primary_particle);

private:

  art::InputTag _strawDigiMCTag;
  art::InputTag _crvDigiMCTag;
  std::vector<art::InputTag> _simParticleTags;
  std::vector<art::InputTag> _extraStepPointMCTags;
  std::vector<art::InputTag> _surfaceStepTags;
  art::InputTag _caloClusterMCTag;
  std::vector<art::InputTag> _crvCoincClusterMCTags;
  art::InputTag _primaryParticleTag;
  art::InputTag _mcTrajectoryTag;
  bool _keepAllGenParticles;
  art::InputTag _strawDigiMCIndexMapTag;
  art::InputTag _crvDigiMCIndexMapTag;
  std::vector<art::InputTag> _caloShowerStepTags;
  art::InputTag _caloShowerSimTag;
  art::InputTag _caloShowerROTag;
  bool _rekeySimParticleCollection;
  std::vector<art::InputTag> _crvStepsToKeep;

  // handles to the old collections
  art::Handle<StrawDigiMCCollection> _strawDigiMCsHandle;
  art::Handle<CrvDigiMCCollection> _crvDigiMCsHandle;
  art::Handle<CaloShowerSimCollection> _caloShowerSimsHandle;
  art::Handle<CaloShowerROCollection> _CaloShowerROsHandle;

  // unique_ptrs to the new output collections
  std::unique_ptr<StrawDigiMCCollection> _newStrawDigiMCs;
  std::unique_ptr<CrvDigiMCCollection> _newCrvDigiMCs;
  std::map<InstanceLabel, std::unique_ptr<StepPointMCCollection> > _newStepPointMCs;
  std::unique_ptr<StrawGasStepCollection> _newStrawGasSteps;
  std::unique_ptr<CrvStepCollection> _newCrvSteps;
  std::unique_ptr<SimParticleCollection> _newSimParticles;
  std::unique_ptr<GenParticleCollection> _newGenParticles;
  std::unique_ptr<CaloShowerStepCollection> _newCaloShowerSteps;
  std::unique_ptr<CaloShowerSimCollection> _newCaloShowerSims;
  std::unique_ptr<CaloShowerROCollection> _newCaloShowerROs;
  std::unique_ptr<SurfaceStepCollection> _newSurfaceSteps;

  // for StepPointMCs, SimParticles and GenParticles we also need reference their new locations with art::Ptrs and so need their ProductIDs and Getters
  std::map<InstanceLabel, art::ProductID> _newStepPointMCsPID;
  std::map<InstanceLabel, const art::EDProductGetter*> _newStepPointMCGetter;
  art::ProductID _newStrawGasStepsPID;
  const art::EDProductGetter* _newStrawGasStepGetter;
  art::ProductID _newCrvStepsPID;
  const art::EDProductGetter* _newCrvStepGetter;
  art::ProductID _newSimParticlesPID;
  const art::EDProductGetter* _newSimParticleGetter;
  art::ProductID _newGenParticlesPID;
  const art::EDProductGetter* _newGenParticleGetter;
  art::ProductID _newCaloShowerStepsPID;
  const art::EDProductGetter* _newCaloShowerStepGetter;
  std::map<art::ProductID, const art::EDProductGetter*> _oldCaloShowerStepGetter;
  art::ProductID _newCaloHitMCsPID;
  const art::EDProductGetter* _newCaloHitMCGetter;
  art::ProductID _newSurfaceStepsPID;
  const art::EDProductGetter* _newSurfaceStepGetter;

  // record the SimParticles that we are keeping so we can use compressSimParticleCollection to do all the work for us
  std::map<art::ProductID, SimParticleSet> _simParticlesToKeep;

  std::vector<InstanceLabel> _newStepPointMCInstances;

  // Optional parameters for reco output
  mu2e::IndexMap _strawDigiMCIndexMap;
  mu2e::IndexMap _crvDigiMCIndexMap;
  art::Handle<CaloClusterMCCollection> _caloClusterMCsHandle;
  std::unique_ptr<CaloClusterMCCollection> _newCaloClusterMCs;
  std::unique_ptr<CaloHitMCCollection> _newCaloHitMCs;
  std::vector<art::Handle<CrvCoincidenceClusterMCCollection>> _crvCoincClusterMCsHandles;
  std::vector<std::unique_ptr<CrvCoincidenceClusterMCCollection>> _newCrvCoincClusterMCs;
  art::Handle<PrimaryParticle> _primaryParticleHandle;
  std::unique_ptr<PrimaryParticle> _newPrimaryParticle;

  // other optional parameters
  art::Handle<MCTrajectoryCollection> _mcTrajectoriesHandle;
  std::unique_ptr<MCTrajectoryCollection> _newMCTrajectories;

  // For CrvDigiMCs, there's a chance that the same StepPointMC will go into multiple CrvDigiMCs
  // This module didn't take this into account initially and so the same StepPointMC was being written out multiple times
  // This std::set is used to make sure that this doesn't happen
  std::set<art::Ptr<CrvStep> > _crvStepsSeen;
  std::map<art::Ptr<CrvStep>, art::Ptr<CrvStep> > _crvStepsMap;

  bool _noCompression;

  // if the map::at fails, produce a useful error message
  inline art::Ptr<SimParticle>& safeRemapRef(SimParticleRemapping& remap, art::Ptr<SimParticle> const& key, int line) const {
    auto it = remap.find(key);
    if(it == remap.end()) {
      throw cet::exception("CompressDigiMCs::safeRemapRef")
        << "remap key "<< key.id() <<" not found at line " << line << "\n";
    }
    return it->second;
  }

};


mu2e::CompressDigiMCs::CompressDigiMCs(const Parameters& conf)
  : art::EDProducer(conf),
    _strawDigiMCTag(conf().strawDigiMCTag()),
    _crvDigiMCTag(conf().crvDigiMCTag()),
    _simParticleTags(conf().simParticleTags()),
    _extraStepPointMCTags(conf().extraStepPointMCTags()),
    _surfaceStepTags(conf().surfaceStepTags()),
    _caloClusterMCTag(conf().caloClusterMCTag()),
    _crvCoincClusterMCTags(conf().crvCoincClusterMCTags()),
    _primaryParticleTag(conf().primaryParticleTag()),
    _mcTrajectoryTag(conf().mcTrajectoryTag()),
    _keepAllGenParticles(conf().keepAllGenParticles()),
  _strawDigiMCIndexMapTag(conf().strawDigiMCIndexMapTag()),
  _crvDigiMCIndexMapTag(conf().crvDigiMCIndexMapTag()),
  _caloShowerStepTags(conf().caloShowerStepTags()),
  _caloShowerSimTag(conf().caloShowerSimTag()),
  _caloShowerROTag(conf().caloShowerROTag()),
  _rekeySimParticleCollection(conf().rekeySimParticleCollection()),
    _crvStepsToKeep(conf().crvStepsToKeep()),
    _noCompression(conf().noCompression())
{
  // Call appropriate produces<>() functions here.
  produces<StrawDigiMCCollection>();
  produces<CrvDigiMCCollection>();

  for (std::vector<art::InputTag>::const_iterator i_tag = _extraStepPointMCTags.begin(); i_tag != _extraStepPointMCTags.end(); ++i_tag) {
    _newStepPointMCInstances.push_back( (*i_tag).instance() );
  }

  for (const auto& i_instance : _newStepPointMCInstances) {
    produces<StepPointMCCollection>( i_instance );
  }
  produces<StrawGasStepCollection>();
  produces<CrvStepCollection>();
  produces<SimParticleCollection>();
  produces<GenParticleCollection>();
  produces<SurfaceStepCollection>();

  // Two possible compressions for calorimeter
  if (_caloShowerStepTags.size() != 0) {
    produces<CaloShowerStepCollection>();
    produces<CaloShowerSimCollection>();
    produces<CaloShowerROCollection>();
  }
  if (_caloClusterMCTag != "") {
    produces<CaloClusterMCCollection>();
    produces<CaloHitMCCollection>();
  }
  for (const auto& crvCoincClusterMCTag : _crvCoincClusterMCTags) {
    produces<CrvCoincidenceClusterMCCollection>(crvCoincClusterMCTag.label());
    _crvCoincClusterMCsHandles.push_back(art::Handle<CrvCoincidenceClusterMCCollection>());
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
  _newStrawGasSteps = std::unique_ptr<StrawGasStepCollection>(new StrawGasStepCollection);
  _newStrawGasStepsPID = event.getProductID<StrawGasStepCollection>();
  _newStrawGasStepGetter = event.productGetter(_newStrawGasStepsPID);

  _newCrvSteps = std::unique_ptr<CrvStepCollection>(new CrvStepCollection);
  _newCrvStepsPID = event.getProductID<CrvStepCollection>();
  _newCrvStepGetter = event.productGetter(_newCrvStepsPID);

  _newSurfaceSteps = std::unique_ptr<SurfaceStepCollection>(new SurfaceStepCollection);
  _newSurfaceStepsPID = event.getProductID<SurfaceStepCollection>();
  _newSurfaceStepGetter = event.productGetter(_newSurfaceStepsPID);

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

    if (_keepAllGenParticles || _noCompression) {
      // Add all the SimParticles that are also GenParticles
      for (const auto& i_oldSimParticle : *oldSimParticles) {
        const cet::map_vector_key& key = i_oldSimParticle.first;
        const SimParticle& i_oldSim = i_oldSimParticle.second;
        if (i_oldSim.genParticle().isNonnull()) {
          keepSimParticle(art::Ptr<SimParticle>(i_product_id, key.asUint(), i_product_getter));
          ++n_gen_particles_to_keep;
        }
        if (_noCompression) { // while we're going through the SimParticleCollection, just add everything if we want no compression
          keepSimParticle(art::Ptr<SimParticle>(i_product_id, key.asUint(), i_product_getter));
        }
      }
    }
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
    _newCaloHitMCs = std::unique_ptr<CaloHitMCCollection>(new CaloHitMCCollection);
    _newCaloHitMCsPID = event.getProductID<CaloHitMCCollection>();
    _newCaloHitMCGetter = event.productGetter(_newCaloHitMCsPID);
  }
  // If we have a CrvCoincClusterMC collection, use that
  if (_crvCoincClusterMCTags.size() > 0) {
    _newCrvCoincClusterMCs.clear();
    for (size_t i_tag = 0; i_tag < _crvCoincClusterMCTags.size(); ++i_tag) {
      const auto & crvCoincClusterMCTag = _crvCoincClusterMCTags.at(i_tag);
      event.getByLabel(crvCoincClusterMCTag, _crvCoincClusterMCsHandles.at(i_tag));
      _newCrvCoincClusterMCs.push_back(std::unique_ptr<CrvCoincidenceClusterMCCollection>(new CrvCoincidenceClusterMCCollection));
    }
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
    if (_strawDigiMCIndexMapTag == "" || in_index_map || _noCompression) {
      copyStrawDigiMC(i_strawDigiMC);
    }
  }

  // Only check for this if we are not reducing the number of StrawDigiMCs
  if ((_strawDigiMCIndexMapTag == "" || _noCompression) && strawDigiMCs.size() != _newStrawDigiMCs->size()) {
    throw cet::exception("CompressDigiMCs") << "The number of StrawDigiMCs before and after compression does not match ("
                                            << strawDigiMCs.size() << " != " << _newStrawDigiMCs->size() << ")" << std::endl;
  }


  if (_crvDigiMCTag != "") {
    _crvStepsSeen.clear();
    _crvStepsMap.clear();

    event.getByLabel(_crvDigiMCTag, _crvDigiMCsHandle);
    const auto& crvDigiMCs = *_crvDigiMCsHandle;
    for (size_t i = 0; i < crvDigiMCs.size(); ++i) {
      const auto& i_crvDigiMC = crvDigiMCs.at(i);
      mu2e::FullIndex full_i = i;
      bool in_index_map = false;
      if (_crvDigiMCIndexMapTag != "") {
        in_index_map = _crvDigiMCIndexMap.checkInMap(full_i);
      }
      if (_crvDigiMCIndexMapTag == "" || in_index_map || _noCompression) {
        copyCrvDigiMC(i_crvDigiMC);
      }
    }
    // Only check for this if we are not reducing the number of CrvDigiMCs
    if ((_crvDigiMCIndexMapTag == "" || _noCompression) && crvDigiMCs.size() != _newCrvDigiMCs->size()) {
      throw cet::exception("CompressDigiMCs") << "The number of CrvDigiMCs before and after compression does not match ("
                                              << crvDigiMCs.size() << " != " << _newCrvDigiMCs->size() << ")" << std::endl;
    }

    // Sometimes we want to keep all CrvSteps regardless of whether they are in a CrvDigiMCs.
    // Here we loop through and add any that were not already included
    for (const auto& crvStepsTag : _crvStepsToKeep) {
      const auto& oldCrvStepsHandle = event.getValidHandle<mu2e::CrvStepCollection>(crvStepsTag);

      // Loop through the CrvSteps we want to keep
      //
      // So far, we have been keeping track of which CrvSteps we have seen through the Ptrs that point to them
      //
      // Since these CrvSteps might not be being referred to by any object, we need to create a fake ptr to make sure we haven't seen it before
      // so let's get the product id and getter so we can construct it
      art::ProductID old_crv_step_product_id = oldCrvStepsHandle.id();
      const art::EDProductGetter* old_crv_step_product_getter = event.productGetter(old_crv_step_product_id);

      for (CrvStepCollection::const_iterator i_crvStep = oldCrvStepsHandle->begin(); i_crvStep != oldCrvStepsHandle->end(); ++i_crvStep) {
        const auto& crvStep = *i_crvStep; // convert from iterator to actual object

        const auto& fake_old_ptr = art::Ptr<CrvStep>(old_crv_step_product_id, i_crvStep - oldCrvStepsHandle->begin(), old_crv_step_product_getter);
        if (_crvStepsSeen.insert(fake_old_ptr).second == true) { // if we have inserted this CrvStepPtrs (i.e. it hasn't already been seen)
          art::Ptr<CrvStep> newStepPtr = copyCrvStep(crvStep);
          _crvStepsMap[fake_old_ptr] = newStepPtr; // need to keep track of these
        }
      }
    }
  }

  // Two possible compressions for calorimeter
  // The first just takes the CaloShowerSteps, CaloShowerSims and CaloShowerROs and reassigns Ptrs (i.e. no actual compression....)
  if (_caloShowerStepTags.size() != 0) {
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

    _newCaloShowerROs = std::unique_ptr<CaloShowerROCollection>(new CaloShowerROCollection);
    event.getByLabel(_caloShowerROTag, _CaloShowerROsHandle);
    const auto& CaloShowerROs = *_CaloShowerROsHandle;
    for (const auto& i_CaloShowerRO : CaloShowerROs) {
      copyCaloShowerRO(i_CaloShowerRO, caloShowerStepRemap);
    }
  }

  // The second uses CaloClusterMCs and only keeps SimParticles that have been assigned to those
  if (_caloClusterMCTag != "") {
    const auto& caloClusterMCs = *_caloClusterMCsHandle;
    for (const auto& i_caloClusterMC : caloClusterMCs) {
      copyCaloClusterMC(i_caloClusterMC);
    }
  }

  // Optional CrvCoincidenceClusterMCs
  for (size_t i_tag = 0; i_tag < _crvCoincClusterMCTags.size(); ++i_tag) {
    const auto& crvCoincClusterMCs = *_crvCoincClusterMCsHandles.at(i_tag);
    for (const auto& i_crvCoincClusterMC : crvCoincClusterMCs) {
      copyCrvCoincClusterMC(i_crvCoincClusterMC, i_tag);
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
        if (!_noCompression) { // if we want to compress
          const SimParticleSet& alreadyKeptSimParts = simPartsToKeep.second;
          for (const auto& alreadyKeptSimPart : alreadyKeptSimParts) {
            if (stepPointMC.simParticle() == alreadyKeptSimPart) {
              copyStepPointMC(stepPointMC, (*i_tag).instance() );
            }
          }
        }
        else { // if we don't want to compress
          copyStepPointMC(stepPointMC, (*i_tag).instance() );
        }
      }
    }
  }

  // Now compress the SurfaceSteps
  for (const auto& surfaceStepsTag : _surfaceStepTags) {

    const auto& oldSurfaceStepsHandle = event.getValidHandle<mu2e::SurfaceStepCollection>(surfaceStepsTag);
    for (const auto& surfaceStep : *oldSurfaceStepsHandle) {
      for (const auto& simPartsToKeep : _simParticlesToKeep) {
        const art::ProductID& oldProdID = simPartsToKeep.first;
        if (surfaceStep.simParticle().id() != oldProdID) {
          continue;
        }
        if (!_noCompression) { // if we want to compress
          const SimParticleSet& alreadyKeptSimParts = simPartsToKeep.second;
          for (const auto& alreadyKeptSimPart : alreadyKeptSimParts) {
            if (surfaceStep.simParticle() == alreadyKeptSimPart) {
              copySurfaceStep(surfaceStep);
            }
          }
        }
        else { // if we don't want to compress
          copySurfaceStep(surfaceStep);
        }
      }
    }
  }
  // Now compress the SimParticleCollections into their new collections
  KeyRemap keyRemap;
  SimParticleRemapping remap;
  unsigned int keep_size = 0;
  for (std::vector<art::InputTag>::const_iterator i_tag = _simParticleTags.begin(); i_tag != _simParticleTags.end(); ++i_tag) {
    keyRemap.clear();
    const auto& oldSimParticles = event.getValidHandle<SimParticleCollection>(*i_tag);
    art::ProductID i_product_id = oldSimParticles.id();
    SimParticleSelector simPartSelector(_simParticlesToKeep[i_product_id]);
    keep_size += _simParticlesToKeep[i_product_id].size();
    if (_rekeySimParticleCollection) {
      compressSimParticleCollection(_newSimParticlesPID, _newSimParticleGetter, *oldSimParticles,
                                    simPartSelector, *_newSimParticles, &keyRemap);
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
        auto it = keyRemap.find(oldKey);
        if(it == keyRemap.end()) {
          throw cet::exception("CompressDigiMCs::badKeyRemap")
            << "keyRemap key "<< oldKey <<" not found\n";
        }
        newKey = it->second;
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
  if ((_keepAllGenParticles || _noCompression) && _newGenParticles->size() != n_gen_particles_to_keep) {
    throw cet::exception("CompressDigiMCs") << "Number of GenParticles in output collection does not match the number of GenParticles we wanted to keep (" << n_gen_particles_to_keep << " != " << _newGenParticles->size() << ")" << std::endl;
  }


  // Now update all objects with SimParticlePtrs

   // Update the StepPointMCs
  for (const auto& i_instance : _newStepPointMCInstances) {
    for (auto& i_stepPointMC : *_newStepPointMCs.at(i_instance)) {
      art::Ptr<SimParticle> newSimPtr = safeRemapRef(remap,i_stepPointMC.simParticle(),__LINE__);
      i_stepPointMC.simParticle() = newSimPtr;
    }
  }

  // Update SurfaceSteps
  for (auto& i_surfaceStep : *_newSurfaceSteps) {
    art::Ptr<SimParticle> newSimPtr = safeRemapRef(remap,i_surfaceStep.simParticle(),__LINE__);
    i_surfaceStep.simParticle() = newSimPtr;
  }

  // Update the StrawGasSteps
  for (auto& i_strawGasStep : *_newStrawGasSteps) {
    art::Ptr<SimParticle> newSimPtr = safeRemapRef(remap,i_strawGasStep.simParticle(),__LINE__);
    i_strawGasStep.simParticle() = newSimPtr;
  }

  // Update the CrvSteps
  if (_crvDigiMCTag != "") {
    for (auto& i_crvStep : *_newCrvSteps) {
      art::Ptr<SimParticle> newSimPtr = safeRemapRef(remap,i_crvStep.simParticle(),__LINE__);
      i_crvStep.simParticle() = newSimPtr;
    }
  }

  if (_caloShowerStepTags.size() != 0) {
    // Update the CaloShowerSteps
    for (auto& i_caloShowerStep : *_newCaloShowerSteps) {
      art::Ptr<SimParticle> newSimPtr = safeRemapRef(remap,i_caloShowerStep.simParticle(),__LINE__);
      i_caloShowerStep.setSimParticle(newSimPtr);
    }
  }

  // NB: copied CaloShowerSims are broken as they refer to the original Step

  if (_caloClusterMCTag != "") {
    for (auto& i_caloHitMC : *_newCaloHitMCs) {
      for (auto& i_caloMCEDep : i_caloHitMC.energyDeposits()) {
        i_caloMCEDep.resetSim(safeRemapRef(remap,i_caloMCEDep.sim(),__LINE__));
      }
    }
  }

  // Update the CrvDigiMCs
  for (auto& i_crvDigiMC : *_newCrvDigiMCs) {
    art::Ptr<SimParticle> oldSimPtr = i_crvDigiMC.GetSimParticle();
    art::Ptr<SimParticle> newSimPtr;
    if (oldSimPtr.isNonnull()) { // if the old CrvDigiMC doesn't have a null ptr for the SimParticle...
      newSimPtr = safeRemapRef(remap,oldSimPtr,__LINE__);
    }
    else {
      newSimPtr = art::Ptr<SimParticle>();
    }
    i_crvDigiMC.setSimParticle(newSimPtr);
  }
  // Update CrvCoincClusterMCs if needs be
  for (size_t i_tag = 0; i_tag < _crvCoincClusterMCTags.size(); ++i_tag) {
    for (auto& i_crvCoincClusterMC : *_newCrvCoincClusterMCs.at(i_tag)) {
      if (i_crvCoincClusterMC.HasMCInfo()) {
        for (auto& i_pulseInfo : i_crvCoincClusterMC.GetModifiablePulses()) {
          art::Ptr<SimParticle> oldSimPtr = i_pulseInfo._simParticle;
          if (oldSimPtr.isNonnull()) { // sometimes CrvCoincedenceClusterMC has DigiMC but no CrvStep or SimParticle (i.e. dark counts)
            art::Ptr<SimParticle> newSimPtr = safeRemapRef(remap,oldSimPtr,__LINE__);
            i_pulseInfo._simParticle = newSimPtr;
          }
        }

        art::Ptr<SimParticle> oldSimPtr = i_crvCoincClusterMC.GetMostLikelySimParticle();
        art::Ptr<SimParticle> newSimPtr = safeRemapRef(remap,oldSimPtr,__LINE__);
        i_crvCoincClusterMC.SetMostLikelySimParticle(newSimPtr);
      }
    }
  }
  // Update PrimaryParticle if needs be
  if (_primaryParticleTag != "") {
    for (auto& i_simPartPtr : _newPrimaryParticle->modifySimParticles()) {
      i_simPartPtr = safeRemapRef(remap,i_simPartPtr,__LINE__);
    }
  }
  // Create new MC Trajectory collection
  if (_mcTrajectoryTag != "") {
    for (const auto& i_mcTrajectory : *_mcTrajectoriesHandle) {
      art::Ptr<SimParticle> oldSimPtr = i_mcTrajectory.first;
      if (remap.find(oldSimPtr) != remap.end()) {
        _newMCTrajectories->insert(std::pair<art::Ptr<SimParticle>, mu2e::MCTrajectory>(safeRemapRef(remap,oldSimPtr,__LINE__), i_mcTrajectory.second));
      }
    }
  }

  // Now add everything to the event
  for (const auto& i_instance : _newStepPointMCInstances) {
    event.put(std::move(_newStepPointMCs.at(i_instance)), i_instance);
  }
  event.put(std::move(_newSurfaceSteps));
  event.put(std::move(_newStrawDigiMCs));
  event.put(std::move(_newStrawGasSteps));
  event.put(std::move(_newCrvSteps));
  event.put(std::move(_newCrvDigiMCs));

  event.put(std::move(_newSimParticles));
  event.put(std::move(_newGenParticles));

  if (_caloShowerStepTags.size() != 0) {
    event.put(std::move(_newCaloShowerSteps));
    event.put(std::move(_newCaloShowerSims));
    event.put(std::move(_newCaloShowerROs));
  }

  if (_caloClusterMCTag != "") {
    event.put(std::move(_newCaloClusterMCs));
    event.put(std::move(_newCaloHitMCs));
  }

  for (size_t i_tag = 0; i_tag < _crvCoincClusterMCTags.size(); ++i_tag) {
    const auto& crvCoincClusterMCTag = _crvCoincClusterMCTags.at(i_tag);
    event.put(std::move(_newCrvCoincClusterMCs.at(i_tag)), crvCoincClusterMCTag.label());
  }
  if (_primaryParticleTag != "") {
    event.put(std::move(_newPrimaryParticle));
  }
  if (_mcTrajectoryTag != "") {
    event.put(std::move(_newMCTrajectories));
  }
}

void mu2e::CompressDigiMCs::copyStrawDigiMC(const mu2e::StrawDigiMC& old_straw_digi_mc) {

  StrawGasStepRemap step_remap;

  // Need to update the Ptrs for the StepPointMCs
  StrawDigiMC::SGSPA newTriggerStepPtr;
  for(int i_end=0;i_end<StrawEnd::nends;++i_end){
    StrawEnd::End end = static_cast<StrawEnd::End>(i_end);

    const auto& old_step_point = old_straw_digi_mc.strawGasStep(end);
    const auto& newStepPtrIter = step_remap.find(old_step_point);
    if (newStepPtrIter == step_remap.end()) {
      if (old_step_point.isAvailable()) {
        step_remap[old_step_point] = copyStrawGasStep( *old_step_point);
      }
      else { // this is a null Ptr but it should be added anyway to keep consistency (not expected for StrawDigis)
        step_remap[old_step_point] = old_step_point;
      }
    }
    art::Ptr<StrawGasStep> new_step_point = step_remap.at(old_step_point);
    newTriggerStepPtr[i_end] = new_step_point;
  }
  StrawDigiMC new_straw_digi_mc(old_straw_digi_mc, newTriggerStepPtr); // copy everything except the Ptrs from the old StrawDigiMC
  _newStrawDigiMCs->push_back(new_straw_digi_mc);
}

void mu2e::CompressDigiMCs::copyCrvDigiMC(const mu2e::CrvDigiMC& old_crv_digi_mc) {

  // Need to update the Ptrs for the StepPointMCs
  std::vector<art::Ptr<CrvStep> > newStepPtrs;
  for (const auto& i_step_mc : old_crv_digi_mc.GetCrvSteps()) {
    if (i_step_mc.isAvailable()) {
      if (_crvStepsSeen.insert(i_step_mc).second == true) { // if we have inserted this CrvStepPtrs (i.e. it hasn't already been seen)
        art::Ptr<CrvStep> newStepPtr = copyCrvStep(*i_step_mc);
        newStepPtrs.push_back(newStepPtr);
        _crvStepsMap[i_step_mc] = newStepPtr;
      }
      else {
        auto it = _crvStepsMap.find(i_step_mc);
        if(it == _crvStepsMap.end()) {
          throw cet::exception("CompressDigiMCs::copyCrvDigiMC")
            << "remap key "<< i_step_mc.id() <<" not found\n";
        }
        newStepPtrs.push_back(it->second);
      }
    }
    else { // this is a null Ptr but it should be added anyway to keep consistency (expected for CrvDigis)
      newStepPtrs.push_back(i_step_mc);
    }
  }

  CrvDigiMC new_crv_digi_mc(old_crv_digi_mc);
  new_crv_digi_mc.setCrvSteps(newStepPtrs);

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
    auto it = remap.find(i_caloShowerStepPtr);
    if(it==remap.end()) {
      throw cet::exception("CompressDigiMCs::copyCaloShowerSim")
        << "remap key "<< i_caloShowerStepPtr.id() <<" not found\n";
    }
    newCaloShowerStepPtrs.push_back(it->second);
  }

  CaloShowerSim new_calo_shower_sim = old_calo_shower_sim;
  new_calo_shower_sim.setCaloShowerSteps(newCaloShowerStepPtrs);

  _newCaloShowerSims->push_back(new_calo_shower_sim);
}

void mu2e::CompressDigiMCs::copyCaloShowerRO(const mu2e::CaloShowerRO& old_calo_shower_step_ro, const CaloShowerStepRemap& remap) {

  const auto& caloShowerStepPtr = old_calo_shower_step_ro.caloShowerStep();
  CaloShowerRO new_calo_shower_step_ro = old_calo_shower_step_ro;
  auto it = remap.find(caloShowerStepPtr);
  if(it == remap.end()) {
    throw cet::exception("CompressDigiMCs::copyCaloShowerRO")
      << "remap key "<< caloShowerStepPtr.id() <<" not found\n";
  }
  new_calo_shower_step_ro.setCaloShowerStep(it->second);

  _newCaloShowerROs->push_back(new_calo_shower_step_ro);
}

art::Ptr<mu2e::CaloHitMC> mu2e::CompressDigiMCs::copyCaloHitMC(const mu2e::CaloHitMC& old_calo_hit_mc) {
  for (const auto& caloEDepMC : old_calo_hit_mc.energyDeposits()) {
    keepSimParticle(caloEDepMC.sim());
  }
  CaloHitMC new_calo_hit_mc(old_calo_hit_mc);
  _newCaloHitMCs->push_back(new_calo_hit_mc);
  return art::Ptr<CaloHitMC>(_newCaloHitMCsPID, _newCaloHitMCs->size()-1, _newCaloHitMCGetter);
}

void mu2e::CompressDigiMCs::copyCaloClusterMC(const mu2e::CaloClusterMC& old_calo_cluster_mc) {
// first, remake the hits
  std::vector<art::Ptr<CaloHitMC>> newhitptrs;
  for(auto const& oldptr : old_calo_cluster_mc.caloHitMCs()){
    newhitptrs.push_back(copyCaloHitMC(*oldptr));
  }
// then, rebuild the cluster
  CaloClusterMC new_calo_cluster_mc(newhitptrs);
  _newCaloClusterMCs->push_back(new_calo_cluster_mc);
}

void mu2e::CompressDigiMCs::copyCrvCoincClusterMC(const mu2e::CrvCoincidenceClusterMC& old_crv_coinc_cluster_mc, size_t i_tag) {

  if (old_crv_coinc_cluster_mc.HasMCInfo()) { // sometimes CrvCoincidenceClusterMC doesn't have MC information e.g. when it is a noise hit
    for (const auto& i_pulseInfo : old_crv_coinc_cluster_mc.GetPulses()) {
      if (i_pulseInfo._simParticle.isNonnull()) { // sometimes CrvCoincedenceClusterMC has DigiMC but no CrvStep or SimParticle (i.e. dark counts)
        keepSimParticle(i_pulseInfo._simParticle);
      }
    }
    keepSimParticle(old_crv_coinc_cluster_mc.GetMostLikelySimParticle());
  }

  CrvCoincidenceClusterMC new_crv_coinc_cluster_mc(old_crv_coinc_cluster_mc);
  _newCrvCoincClusterMCs.at(i_tag)->push_back(new_crv_coinc_cluster_mc);
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

art::Ptr<mu2e::SurfaceStep> mu2e::CompressDigiMCs::copySurfaceStep(const mu2e::SurfaceStep& old_step) {

  keepSimParticle(old_step.simParticle());

  SurfaceStep new_step(old_step);
  _newSurfaceSteps->push_back(new_step);

  return art::Ptr<SurfaceStep>(_newSurfaceStepsPID, _newSurfaceSteps->size()-1, _newSurfaceStepGetter);
}

art::Ptr<mu2e::StrawGasStep> mu2e::CompressDigiMCs::copyStrawGasStep(const mu2e::StrawGasStep& old_step) {

  keepSimParticle(old_step.simParticle());

  StrawGasStep new_step(old_step);
  _newStrawGasSteps->push_back(new_step);

  return art::Ptr<StrawGasStep>(_newStrawGasStepsPID, _newStrawGasSteps->size()-1, _newStrawGasStepGetter);
}

art::Ptr<mu2e::CrvStep> mu2e::CompressDigiMCs::copyCrvStep(const mu2e::CrvStep& old_step) {

  keepSimParticle(old_step.simParticle());

  CrvStep new_step(old_step);
  _newCrvSteps->push_back(new_step);

  return art::Ptr<CrvStep>(_newCrvStepsPID, _newCrvSteps->size()-1, _newCrvStepGetter);
}

void mu2e::CompressDigiMCs::keepSimParticle(const art::Ptr<SimParticle>& sim_ptr) {

  // Also need to add all the parents too
  _simParticlesToKeep[sim_ptr.id()].insert(sim_ptr);
  art::Ptr<SimParticle> childPtr = sim_ptr;
  art::Ptr<SimParticle> parentPtr = childPtr->parent();

  while (parentPtr.isNonnull()) {
    _simParticlesToKeep[sim_ptr.id()].insert(parentPtr);
    childPtr = parentPtr;
    parentPtr = parentPtr->parent();
  }
}


DEFINE_ART_MODULE(mu2e::CompressDigiMCs)
