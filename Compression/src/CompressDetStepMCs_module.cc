//
// Compression module that takes detector steps (StrawGasSteps, CaloShowerSteps, and CrvSteps)
// and compresses out unwanted MC information (StepPointMCs, SimParticles and GenParticles).
// Note that the number of detector steps is NOT reduced.
//
// Allowed compression levels for different data products:
// - StrawGasSteps : noCompression
// - CaloShowerSteps : noCompression
// - CrvSteps : noCompression
// - SurfaceSteps : noCompression
// - SimParticles : noCompression, fullCompression
// - StepPointMCs : noCompression, simParticleCompression
// - MCTrajectories : noCompression, simParticleCompression
//
// There is also the concept of "genealogy" compression, which uses the fhicl parameter
// keepNGenerations and takes an integer corresponding to the
// number of generations back you want to keep. The "oldest" SimParticle remaining
// can be identified as "truncated" with the SimParticle::isTruncated() function
// - Note 1: N = -1 means keep all generations (i.e. no compression)
//
// Dec 2020, Andy Edmonds
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Sequence.h"

#include <memory>

#include "Offline/MCDataProducts/inc/StrawGasStep.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/SurfaceStep.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/Mu2eUtilities/inc/compressSimParticleCollection.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/Compression/inc/CompressionLevel.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/MCDataProducts/inc/MCRelationship.hh"

namespace mu2e {
  class CompressDetStepMCs;

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
}


class mu2e::CompressDetStepMCs : public art::EDProducer {
public:

  struct OptionsConfig {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;

    fhicl::Atom<std::string> strawGasStepCompressionLevel{Name("strawGasStepCompressionLevel"), Comment("Compression level for StrawGasSteps")};
    fhicl::Atom<std::string> caloShowerStepCompressionLevel{Name("caloShowerStepCompressionLevel"), Comment("Compression level for CaloShowerSteps")};
    fhicl::Atom<std::string> crvStepCompressionLevel{Name("crvStepCompressionLevel"), Comment("Compression level for CrvSteps")};
    fhicl::Atom<std::string> surfaceStepCompressionLevel{Name("surfaceStepCompressionLevel"), Comment("Compression level for SurfaceSteps")};
    fhicl::Atom<std::string> simParticleCompressionLevel{Name("simParticleCompressionLevel"), Comment("Compression level for SimParticles")};
    fhicl::Atom<std::string> stepPointMCCompressionLevel{Name("stepPointMCCompressionLevel"), Comment("Compression level for StepPointMCs")};
    fhicl::Atom<int> keepNGenerations{Name("keepNGenerations"), Comment("Number of generations to keep in the genealogy")};
    fhicl::Atom<std::string> mcTrajectoryCompressionLevel{Name("mcTrajectoryCompressionLevel"), Comment("Compression level for MCTrajectories")};
  };

  struct Config {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;

    // fhicl parameters for all of our inputs
    fhicl::Atom<art::InputTag> strawGasStepTag{Name("strawGasStepTag"), Comment("InputTag for the StrawGasSteps")};
    fhicl::Atom<art::InputTag> caloShowerStepTag{Name("caloShowerStepTag"), Comment("InputTag for CaloShowerSteps")};
    fhicl::Atom<art::InputTag> crvStepTag{Name("crvStepTag"), Comment("InputTag for CrvSteps")};
    fhicl::Atom<art::InputTag> surfaceStepTag{Name("surfaceStepTag"), Comment("InputTag for SurfaceSteps")};
    fhicl::Sequence<art::InputTag> stepPointMCTags{Name("stepPointMCTags"), Comment("Sequence of InputTags for StepPointMCCollections (e.g. virtualdetector)")};
    fhicl::Sequence<art::InputTag> simParticleTags{Name("simParticleTags"), Comment("InputTags for the SimParticleCollections")};
    fhicl::Atom<art::InputTag> mcTrajectoryTag{Name("mcTrajectoryTag"), Comment("InputTag for the MCTrajectory")};
    fhicl::Atom<int> debugLevel{Name("debugLevel"), Comment("Debug level (0 = no debug output)")};
    fhicl::Table<OptionsConfig> compressionOptions{Name("compressionOptions"), Comment("Compression options for this module")};
  };
  typedef art::EDProducer::Table<Config> Parameters;

  explicit CompressDetStepMCs(const Parameters& conf);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Required functions.
  void produce(art::Event & event) override;

  // Other functions
  void compressStrawGasSteps(const art::Event& event);
  void updateStrawGasSteps();
  void compressCaloShowerSteps(const art::Event& event);
  void updateCaloShowerSteps();
  void compressCrvSteps(const art::Event& event);
  void updateCrvSteps();
  void compressSurfaceSteps(const art::Event& event);
  void updateSurfaceSteps();
  void compressStepPointMCs(const art::Event& event);
  void updateStepPointMCs();
  void compressSimParticles(const art::Event& event);
  void compressGenParticles();
  void compressMCTrajectories(const art::Event& event);
  void updateMCTrajectories();
  void recordSimParticle(const art::Ptr<mu2e::SimParticle>& sim_ptr);
  void checkCompressionLevels();

private:

  art::InputTag _strawGasStepTag;
  art::InputTag _caloShowerStepTag;
  art::InputTag _crvStepTag;
  art::InputTag _surfaceStepTag;
  std::vector<art::InputTag> _stepPointMCTags;
  std::vector<art::InputTag> _simParticleTags;
  art::InputTag _mcTrajectoryTag;

  mu2e::CompressionLevel _strawGasStepCompressionLevel;
  mu2e::CompressionLevel _caloShowerStepCompressionLevel;
  mu2e::CompressionLevel _crvStepCompressionLevel;
  mu2e::CompressionLevel _surfaceStepCompressionLevel;
  mu2e::CompressionLevel _simParticleCompressionLevel;
  mu2e::CompressionLevel _stepPointMCCompressionLevel;
  int _keepNGenerations;
  mu2e::CompressionLevel _mcTrajectoryCompressionLevel;
  int _debugLevel;

  art::ProductToken<mu2e::StrawGasStepCollection> _strawGasStepToken;
  art::ProductToken<mu2e::CaloShowerStepCollection> _caloShowerStepToken;
  art::ProductToken<mu2e::CrvStepCollection> _crvStepToken;
  art::ProductToken<mu2e::SurfaceStepCollection> _surfaceStepToken;
  art::ProductToken<mu2e::MCTrajectoryCollection> _mcTrajectoryToken;

  // unique_ptrs to the new output collections
  std::unique_ptr<mu2e::StrawGasStepCollection> _newStrawGasSteps;
  std::unique_ptr<CaloShowerStepCollection> _newCaloShowerSteps;
  std::unique_ptr<CrvStepCollection> _newCrvSteps;
  std::unique_ptr<SurfaceStepCollection> _newSurfaceSteps;
  std::map<InstanceLabel, std::unique_ptr<StepPointMCCollection> > _newStepPointMCs;
  std::unique_ptr<MCTrajectoryCollection> _newMCTrajectories;
  // temporary storage for MCTrajectories
  std::vector<MCTrajectory> _newMCTrajs;

  // To create art::Ptrs to the new SimParticles and GenParticles,
  // we need their art::ProductIDs and art::EDProductGetters
  std::unique_ptr<mu2e::SimParticleCollection> _newSimParticles;
  art::ProductID _newSimParticlesPID;
  const art::EDProductGetter* _newSimParticleGetter;

  std::unique_ptr<mu2e::GenParticleCollection> _newGenParticles;
  art::ProductID _newGenParticlesPID;
  const art::EDProductGetter* _newGenParticleGetter;

  // record the SimParticles that we are keeping so we can use compressSimParticleCollection to do all the work for us
  std::map<art::ProductID, mu2e::SimParticleSet> _simParticlesToKeep;
  std::map<art::ProductID, mu2e::SimParticleSet> _simParticlesToTruncate;
  SimParticleRemapping _simPtrRemap;

  // if the map::at fails, produce a useful error message
  inline art::Ptr<SimParticle>& safeRemapRef(art::Ptr<SimParticle> const& key, int line) {
    auto it = _simPtrRemap.find(key);
    if(it == _simPtrRemap.end()) {
      throw cet::exception("CompressDetStepMCs::safeRemapRef")
        << "_simPtrRemap key "<< key.id()<<" not found at line " << line << "\n";
    }
    return it->second;
  }

};


mu2e::CompressDetStepMCs::CompressDetStepMCs(const Parameters& conf)
  : art::EDProducer(conf),
    _strawGasStepTag(conf().strawGasStepTag()),
    _caloShowerStepTag(conf().caloShowerStepTag()),
    _crvStepTag(conf().crvStepTag()),
    _surfaceStepTag(conf().surfaceStepTag()),
    _stepPointMCTags(conf().stepPointMCTags()),
    _simParticleTags(conf().simParticleTags()),
    _mcTrajectoryTag(conf().mcTrajectoryTag()),
    _strawGasStepCompressionLevel(mu2e::CompressionLevel::findByName(conf().compressionOptions().strawGasStepCompressionLevel())),
    _caloShowerStepCompressionLevel(mu2e::CompressionLevel::findByName(conf().compressionOptions().caloShowerStepCompressionLevel())),
    _crvStepCompressionLevel(mu2e::CompressionLevel::findByName(conf().compressionOptions().crvStepCompressionLevel())),
    _surfaceStepCompressionLevel(mu2e::CompressionLevel::findByName(conf().compressionOptions().surfaceStepCompressionLevel())),
    _simParticleCompressionLevel(mu2e::CompressionLevel::findByName(conf().compressionOptions().simParticleCompressionLevel())),
  _stepPointMCCompressionLevel(mu2e::CompressionLevel::findByName(conf().compressionOptions().stepPointMCCompressionLevel())),
  _keepNGenerations(conf().compressionOptions().keepNGenerations()),
  _mcTrajectoryCompressionLevel(mu2e::CompressionLevel::findByName(conf().compressionOptions().mcTrajectoryCompressionLevel())),
  _debugLevel(conf().debugLevel()),
  _strawGasStepToken{mayConsume<mu2e::StrawGasStepCollection>(conf().strawGasStepTag())},
  _caloShowerStepToken{mayConsume<mu2e::CaloShowerStepCollection>(conf().caloShowerStepTag())},
  _crvStepToken{mayConsume<mu2e::CrvStepCollection>(conf().crvStepTag())},
  _surfaceStepToken{mayConsume<mu2e::SurfaceStepCollection>(conf().surfaceStepTag())},
  _mcTrajectoryToken{mayConsume<mu2e::MCTrajectoryCollection>(conf().mcTrajectoryTag())}
{
  // Check that we have valid compression levels for this module
  checkCompressionLevels();

  // Call appropriate produces<>() functions here.
  produces<GenParticleCollection>();
  produces<SimParticleCollection>();
  for (const auto& i_tag : _simParticleTags) {
    consumes<SimParticleCollection>(i_tag);
  }

  produces<StrawGasStepCollection>();
  produces<CaloShowerStepCollection>();
  produces<CrvStepCollection>();
  produces<SurfaceStepCollection>();

  for (const auto& i_tag : _stepPointMCTags) {
    consumes<StepPointMCCollection>(i_tag);
    produces<StepPointMCCollection>( i_tag.instance() );
  }
  produces<MCTrajectoryCollection>();
}

void mu2e::CompressDetStepMCs::produce(art::Event & event)
{
  _newStrawGasSteps = std::unique_ptr<StrawGasStepCollection>(new StrawGasStepCollection);
  _newCaloShowerSteps = std::unique_ptr<CaloShowerStepCollection>(new CaloShowerStepCollection);
  _newCrvSteps = std::unique_ptr<CrvStepCollection>(new CrvStepCollection);
  _newSurfaceSteps = std::unique_ptr<SurfaceStepCollection>(new SurfaceStepCollection);

  for (const auto& i_tag : _stepPointMCTags) {
    _newStepPointMCs[i_tag.instance()] = std::unique_ptr<StepPointMCCollection>(new StepPointMCCollection);
  }

  _newSimParticles = std::unique_ptr<SimParticleCollection>(new SimParticleCollection);
  _newSimParticlesPID = event.getProductID<SimParticleCollection>();
  _newSimParticleGetter = event.productGetter(_newSimParticlesPID);

  _newGenParticles = std::unique_ptr<GenParticleCollection>(new GenParticleCollection);
  _newGenParticlesPID = event.getProductID<GenParticleCollection>();
  _newGenParticleGetter = event.productGetter(_newGenParticlesPID);

  _newMCTrajectories = std::unique_ptr<MCTrajectoryCollection>(new MCTrajectoryCollection);

  _simParticlesToKeep.clear();
  _simParticlesToTruncate.clear();
  _simPtrRemap.clear();

  // Compress detector steps and record which SimParticles we want to keep
  if (_strawGasStepTag != "") { compressStrawGasSteps(event); }
  if (_caloShowerStepTag != "") { compressCaloShowerSteps(event); }
  if (_crvStepTag != "") { compressCrvSteps(event); }
  if (_surfaceStepTag != "") { compressSurfaceSteps(event); }
  if (_stepPointMCCompressionLevel == mu2e::CompressionLevel::kNoCompression) {
    // if we are not compressing StepPointMCs, then
    // we want to make sure we record all their SimParticles
    compressStepPointMCs(event);
  }
  if (_mcTrajectoryTag != "") {
    if (_mcTrajectoryCompressionLevel == mu2e::CompressionLevel::kNoCompression) {
      // if we are not compressing MCTrajectories, then
      // we want to make sure we record all their SimParticles
      compressMCTrajectories(event);
    }
  }

  // Compress the SimParticles and record their new keys
  compressSimParticles(event);

  // Now that we know which SimParticles we are keeping,
  // we will keep the data products that are associated with these SimParticles
  compressGenParticles();
  if (_stepPointMCCompressionLevel == mu2e::CompressionLevel::kSimParticleCompression) {
    // if we are compressing StepPointMCs based on the SimParticles we are keeping,
    // then compressStepPointMCs now
    compressStepPointMCs(event);
  }
  if (_mcTrajectoryTag != "") {
    if (_mcTrajectoryCompressionLevel == mu2e::CompressionLevel::kSimParticleCompression) {
      // if we are compressing MCTrajectories based on the SimParticles we are keeping,
      // then compress MCTrajectories now
      compressMCTrajectories(event);
    }
  }

  // Update all the data products so that their SimParticlePtrs point to the new collection
  if (_strawGasStepTag != "") { updateStrawGasSteps(); }
  if (_caloShowerStepTag != "") { updateCaloShowerSteps(); }
  if (_crvStepTag != "") { updateCrvSteps(); }
  if (_surfaceStepTag != "") { updateSurfaceSteps(); }
  updateStepPointMCs();
  if (_mcTrajectoryTag != "") { updateMCTrajectories(); }

  // Now add everything to the event
  event.put(std::move(_newStrawGasSteps));
  event.put(std::move(_newCaloShowerSteps));
  event.put(std::move(_newCrvSteps));
  event.put(std::move(_newSurfaceSteps));
  for (const auto& i_tag : _stepPointMCTags) {
    event.put(std::move(_newStepPointMCs.at(i_tag.instance())), i_tag.instance());
  }
  event.put(std::move(_newSimParticles));
  event.put(std::move(_newGenParticles));
  event.put(std::move(_newMCTrajectories));
}

void mu2e::CompressDetStepMCs::compressStrawGasSteps(const art::Event& event) {
  const auto& strawGasStepsHandle = event.getValidHandle(_strawGasStepToken);
  const auto& strawGasSteps = *strawGasStepsHandle;
  if(_debugLevel>0 && strawGasSteps.size()>0) {
    std::cout << "Compressing StrawGasSteps from " << _strawGasStepTag << std::endl;
  }
  for (const auto& i_strawGasStep : strawGasSteps) {
    if (_simParticleCompressionLevel == mu2e::CompressionLevel::kFullCompression) {
      recordSimParticle(i_strawGasStep.simParticle());
    }
    StrawGasStep newStrawGasStep(i_strawGasStep);
    _newStrawGasSteps->push_back(newStrawGasStep);
  }
  if(_strawGasStepCompressionLevel == mu2e::CompressionLevel::kNoCompression) {
    if (_newStrawGasSteps->size() != strawGasSteps.size()) {
      throw cet::exception("CompressDetStepMCs") << "Number of StrawGasSteps in output collection (" << _newStrawGasSteps->size() << ") does not match the number of StrawGasSteps in the input collection (" << strawGasSteps.size() << ") even though no compression has been requested (strawGasStepCompressionLevel = \"" << _strawGasStepCompressionLevel.name() << "\")" << std::endl;
    }
  }
}

void mu2e::CompressDetStepMCs::compressCaloShowerSteps(const art::Event& event) {
  const auto& caloShowerStepsHandle = event.getValidHandle(_caloShowerStepToken);
  const auto& caloShowerSteps = *caloShowerStepsHandle;
  if(_debugLevel>0 && caloShowerSteps.size()>0) {
    std::cout << "Compressing CaloShowerSteps from " << _caloShowerStepTag << std::endl;
  }
  for (const auto& i_caloShowerStep : caloShowerSteps) {
    if (_simParticleCompressionLevel == mu2e::CompressionLevel::kFullCompression) {
      recordSimParticle(i_caloShowerStep.simParticle());
    }
    CaloShowerStep newCaloShowerStep(i_caloShowerStep);
    _newCaloShowerSteps->push_back(newCaloShowerStep);
  }
  if(_caloShowerStepCompressionLevel == mu2e::CompressionLevel::kNoCompression) {
    if (_newCaloShowerSteps->size() != caloShowerSteps.size()) {
      throw cet::exception("CompressDetStepMCs") << "Number of CaloShowerSteps in output collection (" << _newCaloShowerSteps->size() << ") does not match the number of CaloShowerSteps in the input collection (" << caloShowerSteps.size() << ") even though no compression has been requested (caloShowerStepCompressionLevel = \"" << _caloShowerStepCompressionLevel.name() << "\")" << std::endl;
    }
  }
}

void mu2e::CompressDetStepMCs::compressCrvSteps(const art::Event& event) {
  const auto& crvStepsHandle = event.getValidHandle(_crvStepToken);
  const auto& crvSteps = *crvStepsHandle;
  if(_debugLevel>0 && crvSteps.size()>0) {
    std::cout << "Compressing CrvSteps from " << _crvStepTag << std::endl;
  }
  for (const auto& i_crvStep : crvSteps) {
    if (_simParticleCompressionLevel == mu2e::CompressionLevel::kFullCompression) {
      recordSimParticle(i_crvStep.simParticle());
    }
    CrvStep newCrvStep(i_crvStep);
    _newCrvSteps->push_back(newCrvStep);
  }
  if(_crvStepCompressionLevel == mu2e::CompressionLevel::kNoCompression) {
    if (_newCrvSteps->size() != crvSteps.size()) {
      throw cet::exception("CompressDetStepMCs") << "Number of CrvSteps in output collection (" << _newCrvSteps->size() << ") does not match the number of CrvSteps in the input collection (" << crvSteps.size() << ") even though no compression has been requested (crvStepCompressionLevel = \"" << _crvStepCompressionLevel.name() << "\")" << std::endl;
    }
  }
}

void mu2e::CompressDetStepMCs::compressSurfaceSteps(const art::Event& event) {
  const auto& surfaceStepsHandle = event.getValidHandle(_surfaceStepToken);
  const auto& surfaceSteps = *surfaceStepsHandle;
  if(_debugLevel>0 && surfaceSteps.size()>0) {
    std::cout << "Compressing SurfaceSteps from " << _surfaceStepTag << std::endl;
  }
  for (const auto& i_surfaceStep : surfaceSteps) {
    if (_simParticleCompressionLevel == mu2e::CompressionLevel::kFullCompression) {
      recordSimParticle(i_surfaceStep.simParticle());
    }
    SurfaceStep newSurfaceStep(i_surfaceStep);
    _newSurfaceSteps->push_back(newSurfaceStep);
  }
  if(_surfaceStepCompressionLevel == mu2e::CompressionLevel::kNoCompression) {
    if (_newSurfaceSteps->size() != surfaceSteps.size()) {
      throw cet::exception("CompressDetStepMCs") << "Number of SurfaceSteps in output collection (" << _newSurfaceSteps->size() << ") does not match the number of SurfaceSteps in the input collection (" << surfaceSteps.size() << ") even though no compression has been requested (surfaceStepCompressionLevel = \"" << _surfaceStepCompressionLevel.name() << "\")" << std::endl;
    }
  }
}

void mu2e::CompressDetStepMCs::updateStrawGasSteps() {
  for (auto& i_strawGasStep : *_newStrawGasSteps) {
    const auto& oldSimPtr = i_strawGasStep.simParticle();
    art::Ptr<mu2e::SimParticle> newSimPtr = safeRemapRef( oldSimPtr, __LINE__);
    if(_debugLevel>0) {
      std::cout << "Updating SimParticlePtr in StrawGasStep from " << oldSimPtr << " to " << newSimPtr << std::endl;
    }

    i_strawGasStep.simParticle() = newSimPtr;
  }
}

void mu2e::CompressDetStepMCs::updateCaloShowerSteps() {
  for (auto& i_caloShowerStep : *_newCaloShowerSteps) {
    const auto& oldSimPtr = i_caloShowerStep.simParticle();
    art::Ptr<mu2e::SimParticle> newSimPtr = safeRemapRef( oldSimPtr, __LINE__);;
    if(_debugLevel>0) {
      std::cout << "Updating SimParticlePtr in CaloShowerStep from " << oldSimPtr << " to " << newSimPtr << std::endl;
    }

    i_caloShowerStep.setSimParticle(newSimPtr);
  }
}

void mu2e::CompressDetStepMCs::updateCrvSteps() {
  for (auto& i_crvStep : *_newCrvSteps) {
    const auto& oldSimPtr = i_crvStep.simParticle();
    art::Ptr<mu2e::SimParticle> newSimPtr = safeRemapRef( oldSimPtr, __LINE__);
    if(_debugLevel>0) {
      std::cout << "Updating SimParticlePtr in CrvStep from " << oldSimPtr << " to " << newSimPtr << std::endl;
    }

    i_crvStep.simParticle() = newSimPtr;
  }
}

void mu2e::CompressDetStepMCs::updateSurfaceSteps() {
  for (auto& i_surfaceStep : *_newSurfaceSteps) {
    const auto& oldSimPtr = i_surfaceStep.simParticle();
    art::Ptr<mu2e::SimParticle> newSimPtr = safeRemapRef( oldSimPtr, __LINE__);
    if(_debugLevel>0) {
      std::cout << "Updating SimParticlePtr in SurfaceStep from " << oldSimPtr << " to " << newSimPtr << std::endl;
    }

    i_surfaceStep.simParticle() = newSimPtr;
  }
}

void mu2e::CompressDetStepMCs::compressSimParticles(const art::Event& event) {
  // Now compress the SimParticleCollections into their new collections
  KeyRemap* keyRemap = new KeyRemap; // if we have multiple SimParticleCollections, we will need to rekey the SimParticles
  bool rekeySimParticleCollection = false;
  if (_simParticleTags.size() > 1) {
    rekeySimParticleCollection = true;
  }
  unsigned int keep_size = 0;
  for (const auto& i_tag : _simParticleTags) {
    keyRemap->clear();
    const auto& oldSimParticles = event.getValidHandle<SimParticleCollection>(i_tag);
    art::ProductID i_product_id = oldSimParticles.id();
    const art::EDProductGetter* i_prod_getter = event.productGetter(i_product_id);
    if (_simParticleCompressionLevel == CompressionLevel::kNoCompression) {
      // add all the SimParticles
      for (const auto& i_simParticle : *oldSimParticles) {
        art::Ptr<SimParticle> oldSimPtr(i_product_id, i_simParticle.first.asUint(), i_prod_getter);
        recordSimParticle(oldSimPtr);
      }
    }

    SimParticleSelector simPartSelector(_simParticlesToKeep[i_product_id]);
    keep_size += _simParticlesToKeep[i_product_id].size();
    if (rekeySimParticleCollection) {
      compressSimParticleCollection(_newSimParticlesPID, _newSimParticleGetter, *oldSimParticles, simPartSelector, *_newSimParticles, keyRemap);
    }
    else {
      compressSimParticleCollection(_newSimParticlesPID, _newSimParticleGetter, *oldSimParticles, simPartSelector, *_newSimParticles);
    }

    // Fill out the SimParticleRemapping
    for (const auto& i_keptSimPart : _simParticlesToKeep[i_product_id]) {
      cet::map_vector_key oldKey = cet::map_vector_key(i_keptSimPart.key());
      cet::map_vector_key newKey = oldKey;
      if (rekeySimParticleCollection) {
        auto it = keyRemap->find(oldKey);
        if(it == keyRemap->end()) {
          throw cet::exception("CompressDetStepMCs::compressSimParticles") << "Failed to find key "
                                        << oldKey.asUint() << " at line "<< __LINE__  << "\n";
        }
        newKey = it->second;
      }
      _simPtrRemap[i_keptSimPart] = art::Ptr<mu2e::SimParticle>(_newSimParticlesPID, newKey.asUint(), _newSimParticleGetter);
      if (_debugLevel>0) {
        std::cout << "Compressing SimParticle " << i_keptSimPart << " --> " << safeRemapRef( i_keptSimPart, __LINE__) << std::endl;
      }
    }
    if (keep_size != _newSimParticles->size()) {
      throw cet::exception("CompressDetStepMCs") << "Number of SimParticles in output collection (" << _newSimParticles->size() << ") does not match the number of SimParticles we wanted to keep (" << keep_size << ")" << std::endl;
    }
    if (_simParticleCompressionLevel == mu2e::CompressionLevel::kNoCompression) {
      if (_newSimParticles->size() != oldSimParticles->size()) {
        throw cet::exception("CompressDetStepMCs") << "Number of SimParticles in output collection (" << _newSimParticles->size() << ") does not match the number of SimParticles in the input collection (" << oldSimParticles->size() << ") even though no compression has been requested (simParticleCompressionLevel = \"" << _simParticleCompressionLevel.name() << "\")" << std::endl;
      }
    }

    // If we asked for the genealogy to be compressed, we will now end up with some missing links which we need to fix
    // (these should all be within a single input SimParticleCollection)
    if (_keepNGenerations >= 0) {
      // Go through the particles we are keeping and see if any parents are not there
      for (auto& i_keptSimPart : _simParticlesToKeep[i_product_id]) {

        art::Ptr<mu2e::SimParticle> i_childPtr = i_keptSimPart;
        art::Ptr<mu2e::SimParticle> i_parentPtr = i_childPtr->parent();
        while (i_parentPtr) {
          // if the parent will not be in the output collection
          if (_simPtrRemap.find(i_parentPtr) == _simPtrRemap.end()) {
            if (_debugLevel>0) {
              std::cout << "SimParticle " << i_parentPtr << " will not be in output collection because it has been compressed away by genealogy compression" << std::endl;
            }

            _simParticlesToTruncate[i_childPtr.id()].insert(i_childPtr);
            break; // don't go further up the genealogy tree otherwise we will be adding particles
          }
          else {
            // this parent is in the output collection so
            if (_debugLevel>0) {
              std::cout << "SimParticle " << i_parentPtr << " is in the output collection as " << safeRemapRef( i_parentPtr, __LINE__) << std::endl;
            }
          }
          i_childPtr = i_parentPtr;
          i_parentPtr = i_parentPtr->parent();
        }
      }

      // Go through the truncated SimParticles and fix the parent/child links
      for (const auto& i_truncatedSimPart : _simParticlesToTruncate[i_product_id]) {
        //    for (auto& i_simParticle : *_newSimParticles) {
        mu2e::SimParticle& newsim = (*_newSimParticles)[i_truncatedSimPart->id()];//_newSimParticles->at(i_truncatedSimPart.second);//i_simParticle.second;
        // go up genealogy to get the next ancestor that is in the output
        art::Ptr<mu2e::SimParticle> i_ancestorPtr = newsim.parent();
        if (_debugLevel>0) {
          std::cout << "Look for a new parent for particle id " << newsim.id() << " (current parent = " << i_ancestorPtr << ")" << std::endl;
        }
        while (i_ancestorPtr) {
          const auto& findIter = _simPtrRemap.find(i_ancestorPtr);
          if (findIter != _simPtrRemap.end()) {
            newsim.parent() = findIter->second;
            art::Ptr<mu2e::SimParticle> newChildPtr = art::Ptr<mu2e::SimParticle>(_newSimParticlesPID, newsim.id().asUint(), _newSimParticleGetter);
            (*_newSimParticles)[i_ancestorPtr->id()].addDaughter(newChildPtr);
            if (_debugLevel > 0) {
              std::cout << "Because of truncation setting SimParticle (" << newsim.id() << ")'s parent to " << findIter->second << " and adding daughter " << newChildPtr << std::endl;
            }
            break; // don't need to go any further
          }
          else {
            // If we have got to the very first SimParticle (i.e. the one that points to the GenParticle)
            if (i_ancestorPtr->isPrimary()) {
              newsim.genParticle() = i_ancestorPtr->genParticle();// set this particle's GenParticlePtr
              newsim.parent() = art::Ptr<SimParticle>(); // remove the parent
              break; // don't need to go any further
            }
            else { // this is just another step in the genealogy
              i_ancestorPtr = i_ancestorPtr->parent();
            }
          }
        }
      }
    }
  }

  if (_debugLevel > 0) {
    std::cout << "Final SimParticleCollection:" << std::endl;
    for (auto& i_simParticle : *_newSimParticles) {
      mu2e::SimParticle& newsim = i_simParticle.second;
      std::cout << "id = " << i_simParticle.first << ", pdg = " << newsim.pdgId() << ", creation code = " << newsim.creationCode() << ", parent = " << newsim.parent() << ", daughters: ";
      for (const auto& i_daughter : newsim.daughters()) {
        std::cout << i_daughter << " ";
      }
      std::cout << std::endl;
    }
  }
}

void mu2e::CompressDetStepMCs::compressGenParticles() {
  // Loop through the new SimParticles to keep any GenParticles
  if (_debugLevel > 0) {
    std::cout << "Compressing GenParticles..." << std::endl;
  }
  for (auto& i_simParticle : *_newSimParticles) {
    mu2e::SimParticle& newsim = i_simParticle.second;
    if(newsim.genParticle().isNonnull()) { // will crash if not resolvable
      // Copy GenParticle to the new collection
      _newGenParticles->emplace_back(*newsim.genParticle());
      newsim.genParticle() = art::Ptr<mu2e::GenParticle>(_newGenParticlesPID, _newGenParticles->size()-1, _newGenParticleGetter);
      if (_debugLevel > 0) {
        std::cout << "Keeping GenParticle with Ptr " << newsim.genParticle() << std::endl;
      }
    }
  }
}

void mu2e::CompressDetStepMCs::compressStepPointMCs(const art::Event& event) {

  for (const auto& i_tag : _stepPointMCTags) {
    const auto& stepPointMCs = event.getValidHandle<StepPointMCCollection>(i_tag);
    if(_debugLevel>0 && stepPointMCs->size()>0) {
      std::cout << "Compressing StepPointMCs from " << i_tag << std::endl;
    }
    for (const auto& stepPointMC : *stepPointMCs) {
      if (_stepPointMCCompressionLevel == mu2e::CompressionLevel::kSimParticleCompression) {
        for (const auto& simPartsToKeep : _simParticlesToKeep) {
          const art::ProductID& oldProdID = simPartsToKeep.first;
          if (stepPointMC.simParticle().id() != oldProdID) {
            continue;
          }
          const SimParticleSet& alreadyKeptSimParts = simPartsToKeep.second;
          for (const auto& alreadyKeptSimPart : alreadyKeptSimParts) {
            if (stepPointMC.simParticle() == alreadyKeptSimPart) {
              StepPointMC newStepPointMC(stepPointMC);
              _newStepPointMCs.at(i_tag.instance())->push_back(newStepPointMC);
            }
          }
        }
      }
      else if (_stepPointMCCompressionLevel == mu2e::CompressionLevel::kNoCompression) {
        StepPointMC newStepPointMC(stepPointMC);
        _newStepPointMCs.at(i_tag.instance())->push_back(newStepPointMC);
        recordSimParticle(stepPointMC.simParticle());
      }
      else {
        throw cet::exception("CompressDetStepMCs") << "Unrecognized compression level \"" << _stepPointMCCompressionLevel.name() <<"\" for StepPointMCs" << std::endl;
      }
    }
  }
}
void mu2e::CompressDetStepMCs::compressMCTrajectories(const art::Event& event) {

  const auto& mcTrajectories = event.getValidHandle(_mcTrajectoryToken);
  if(_debugLevel>0 && mcTrajectories->size()>0) {
    std::cout << "Compressing MCTrajectories from " << _mcTrajectoryTag << " (size = " << mcTrajectories->size() << ")" << std::endl;
  }
  for (const auto& mcTrajectory : *mcTrajectories) {
    if (_mcTrajectoryCompressionLevel == mu2e::CompressionLevel::kSimParticleCompression) {
      for (const auto& simPartsToKeep : _simParticlesToKeep) {
        const art::ProductID& oldProdID = simPartsToKeep.first;
        if (mcTrajectory.second.sim().id() != oldProdID) {
          continue;
        }
        const SimParticleSet& alreadyKeptSimParts = simPartsToKeep.second;
        for (const auto& alreadyKeptSimPart : alreadyKeptSimParts) {
          if (mcTrajectory.second.sim() == alreadyKeptSimPart) {
            _newMCTrajs.emplace_back(mcTrajectory.second);
            if(_debugLevel>0 ) std::cout << "Inserting new MCTrajectory with key " << mcTrajectory.second.sim() << std::endl;
          }
        }
      }
    }
    else if (_mcTrajectoryCompressionLevel == mu2e::CompressionLevel::kNoCompression) {
      _newMCTrajs.emplace_back(mcTrajectory.second);
      if(_debugLevel>0 ) std::cout << "Inserting new MCTrajectory with key " << mcTrajectory.second.sim() << " and recording it " << std::endl;
      recordSimParticle(mcTrajectory.second.sim());
    }
    else {
      throw cet::exception("CompressDetStepMCs") << "Unrecognized compression level \"" << _mcTrajectoryCompressionLevel.name() <<"\" for MCTrajectories" << std::endl;
    }
  }
}

void mu2e::CompressDetStepMCs::updateStepPointMCs() {
  for (const auto& i_tag : _stepPointMCTags) {
    for (auto& i_stepPointMC : *(_newStepPointMCs.at(i_tag.instance()))) {
      const auto& oldSimPtr = i_stepPointMC.simParticle();
      art::Ptr<mu2e::SimParticle> newSimPtr = safeRemapRef( oldSimPtr, __LINE__);
      if(_debugLevel>0) {
        std::cout << "Updating SimParticlePtr in StepPointMC from " << oldSimPtr << " to " << newSimPtr << std::endl;
      }
      i_stepPointMC.simParticle() = newSimPtr;
    }
  }
}

void mu2e::CompressDetStepMCs::updateMCTrajectories() {
  for (auto& i_mcTrajectory : _newMCTrajs) {
    const auto& oldSimPtr = i_mcTrajectory.sim();
    art::Ptr<mu2e::SimParticle> newSimPtr = safeRemapRef( oldSimPtr, __LINE__);
    if(_debugLevel>0) {
      std::cout << "Updating SimParticlePtr in MCTrajectory from " << oldSimPtr << " to " << newSimPtr << std::endl;
    }
    // overwrite the internal ptr and insert with the new key
    i_mcTrajectory.sim() = newSimPtr;
    _newMCTrajectories->insert(std::make_pair(i_mcTrajectory.sim(), i_mcTrajectory));
  }
  _newMCTrajs.clear(); // clear the temporary storage
}

void mu2e::CompressDetStepMCs::recordSimParticle(const art::Ptr<mu2e::SimParticle>& sim_ptr) {
  // Also need to add all the parents too
  _simParticlesToKeep[sim_ptr.id()].insert(sim_ptr);
  art::Ptr<mu2e::SimParticle> childPtr = sim_ptr;
  art::Ptr<mu2e::SimParticle> parentPtr = childPtr->parent();

  if(_debugLevel>0) {
    std::cout << "Recording SimParticle " << sim_ptr << std::endl;
  }
  while (parentPtr) {
    MCRelationship mcr(sim_ptr, parentPtr);
    if (_keepNGenerations == -1 || ( (mcr.removal() <= _keepNGenerations) && mcr.removal()>=0) ) {
      _simParticlesToKeep[sim_ptr.id()].insert(parentPtr);
      if(_debugLevel>0) {
        std::cout << "and recording its ancestor " << parentPtr << " (NGen = " << (int)mcr.removal() << ")" << std::endl;
      }
    }
    // else if (parentPtr->isPrimary()) { // always keep the very first SimParticle
    //   _simParticlesToKeep[sim_ptr.id()].insert(parentPtr);
    //   if(_debugLevel>0) {
    //     std::cout << "and recording the very first SimParticle " << parentPtr << std::endl;
    //   }
    // }
    else {
      if(_debugLevel>0) {
        std::cout << "and *not* recording its ancestor " << parentPtr << " (NGen = " << (int)mcr.removal() << ")" << std::endl;
      }
    }
    childPtr = parentPtr;
    parentPtr = parentPtr->parent();
  }
}

void mu2e::CompressDetStepMCs::checkCompressionLevels() {
  if (_strawGasStepCompressionLevel != mu2e::CompressionLevel::kNoCompression) {
    throw cet::exception("CompressDetStepMCs") << "This module does not allow StrawGasSteps to be compressed with compression level \"" << _strawGasStepCompressionLevel.name() <<"\"" << std::endl;
  }

  if (_caloShowerStepCompressionLevel != mu2e::CompressionLevel::kNoCompression) {
    throw cet::exception("CompressDetStepMCs") << "This module does not allow CaloShowerSteps to be compressed with compression level \"" << _caloShowerStepCompressionLevel.name() <<"\"" << std::endl;
  }

  if (_crvStepCompressionLevel != mu2e::CompressionLevel::kNoCompression) {
    throw cet::exception("CompressDetStepMCs") << "This module does not allow CrvSteps to be compressed with compression level \"" << _crvStepCompressionLevel.name() <<"\"" << std::endl;
  }

  if (_surfaceStepCompressionLevel != mu2e::CompressionLevel::kNoCompression) {
    throw cet::exception("CompressDetStepMCs") << "This module does not allow SurfaceSteps to be compressed with compression level \"" << _surfaceStepCompressionLevel.name() <<"\"" << std::endl;
  }

  if (_simParticleCompressionLevel != mu2e::CompressionLevel::kNoCompression &&
      _simParticleCompressionLevel != mu2e::CompressionLevel::kFullCompression) {
    throw cet::exception("CompressDetStepMCs") << "This module does not allow SimParticles to be compressed with compression level \"" << _simParticleCompressionLevel.name() <<"\"" << std::endl;
  }

  if (_stepPointMCCompressionLevel != mu2e::CompressionLevel::kNoCompression &&
      _stepPointMCCompressionLevel != mu2e::CompressionLevel::kSimParticleCompression) {
    throw cet::exception("CompressDetStepMCs") << "This module does not allow StepPointMCs to be compressed with compression level \"" << _stepPointMCCompressionLevel.name() <<"\"" << std::endl;
  }

  if (_mcTrajectoryCompressionLevel != mu2e::CompressionLevel::kNoCompression &&
      _mcTrajectoryCompressionLevel != mu2e::CompressionLevel::kSimParticleCompression) {
    throw cet::exception("CompressDetStepMCs") << "This module does not allow MCTrajectories to be compressed with compression level \"" << _mcTrajectoryCompressionLevel.name() <<"\"" << std::endl;
  }
}

DEFINE_ART_MODULE(mu2e::CompressDetStepMCs)
