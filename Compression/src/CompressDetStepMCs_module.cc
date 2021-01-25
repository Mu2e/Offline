//
// Compression module that takes detector steps (StrawGasSteps, CaloShowerSteps, and CrvSteps)
// and compresses out unwanted MC information (StepPointMCs, SimParticles and GenParticles).
// Note that the number of detector steps is NOT reduced.
//
// Allowed compression levels for different data products:
// - StrawGasSteps : noCompression
// - CaloShowerSteps : noCompression
// - CrvSteps : noCompression
// - SimParticles : noCompression, fullCompression
// - StepPointMCs : noCompression, simParticleCompression
// - MCTrajectories : noCompression, simParticleCompression
// - PrimaryParticle : noCompression
//
// There is also the concept of "genealogy" compression, the options are:
// - noCompression : self-explanatory
// - fullCompression : remove all SimParticles between the ones we are keeping and the very first SimParticle
// in all cases the very first SimParticle is kept
//
// Dec 2020, Andy Edmonds
//
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

#include "MCDataProducts/inc/StrawGasStep.hh"
#include "MCDataProducts/inc/CaloShowerStep.hh"
#include "MCDataProducts/inc/CrvStep.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "Mu2eUtilities/inc/compressSimParticleCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Compression/inc/CompressionLevel.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/PrimaryParticle.hh"

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
  struct Config {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;

    // fhicl parameters for all of our inputs
    fhicl::Atom<art::InputTag> strawGasStepTag{Name("strawGasStepTag"), Comment("InputTag for the StrawGasSteps")};
    fhicl::Atom<art::InputTag> caloShowerStepTag{Name("caloShowerStepTag"), Comment("InputTag for CaloShowerSteps")};
    fhicl::Atom<art::InputTag> crvStepTag{Name("crvStepTag"), Comment("InputTag for CrvSteps")};
    fhicl::Sequence<art::InputTag> stepPointMCTags{Name("stepPointMCTags"), Comment("Sequence of InputTags for StepPointMCCollections (e.g. virtualdetector)")};
    fhicl::Atom<art::InputTag> simParticleTag{Name("simParticleTag"), Comment("InputTag for the SimParticleCollection")};
    fhicl::Atom<int> debugLevel{Name("debugLevel"), Comment("Debug level (0 = no debug output)")};
    fhicl::Atom<std::string> strawGasStepCompressionLevel{Name("strawGasStepCompressionLevel"), Comment("Compression level for StrawGasSteps")};
    fhicl::Atom<std::string> caloShowerStepCompressionLevel{Name("caloShowerStepCompressionLevel"), Comment("Compression level for CaloShowerSteps")};
    fhicl::Atom<std::string> crvStepCompressionLevel{Name("crvStepCompressionLevel"), Comment("Compression level for CrvSteps")};
    fhicl::Atom<std::string> simParticleCompressionLevel{Name("simParticleCompressionLevel"), Comment("Compression level for SimParticles")};
    fhicl::Atom<std::string> stepPointMCCompressionLevel{Name("stepPointMCCompressionLevel"), Comment("Compression level for StepPointMCs")};
    fhicl::Atom<std::string> genealogyCompressionLevel{Name("genealogyCompressionLevel"), Comment("Compression level for the genealogy")};
    fhicl::Atom<int> truncatedSimParticleKeyOffset{Name("truncatedSimParticleKeyOffset"), Comment("Offset to use when adding a truncated SimParticle to the SimParticleCollection")};
    fhicl::Atom<art::InputTag> mcTrajectoryTag{Name("mcTrajectoryTag"), Comment("InputTag for the SimParticleCollection")};
    fhicl::Atom<std::string> mcTrajectoryCompressionLevel{Name("mcTrajectoryCompressionLevel"), Comment("Compression level for MCTrajectories")};
    fhicl::Atom<art::InputTag> primaryParticleTag{Name("primaryParticleTag"), Comment("InputTag for the SimParticleCollection")};
    fhicl::Atom<std::string> primaryParticleCompressionLevel{Name("primaryParticleCompressionLevel"), Comment("Compression level for MCTrajectories")};
  };
  typedef art::EDProducer::Table<Config> Parameters;

  explicit CompressDetStepMCs(const Parameters& conf);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CompressDetStepMCs(CompressDetStepMCs const &) = delete;
  CompressDetStepMCs(CompressDetStepMCs &&) = delete;
  CompressDetStepMCs & operator = (CompressDetStepMCs const &) = delete;
  CompressDetStepMCs & operator = (CompressDetStepMCs &&) = delete;

  // Required functions.
  void produce(art::Event & event) override;

  // Other functions
  void compressStrawGasSteps(const art::Event& event);
  void updateStrawGasSteps();
  void compressCaloShowerSteps(const art::Event& event);
  void updateCaloShowerSteps();
  void compressCrvSteps(const art::Event& event);
  void updateCrvSteps();
  void compressStepPointMCs(const art::Event& event);
  void updateStepPointMCs();
  void compressSimParticles(const art::Event& event);
  void compressGenParticles();
  void compressMCTrajectories(const art::Event& event);
  void updateMCTrajectories();
  void compressPrimaryParticle(const art::Event& event);
  void updatePrimaryParticle();
  void recordSimParticle(const art::Ptr<mu2e::SimParticle>& sim_ptr);
  void checkCompressionLevels();

private:

  Config _conf;
  mu2e::CompressionLevel _strawGasStepCompressionLevel;
  mu2e::CompressionLevel _caloShowerStepCompressionLevel;
  mu2e::CompressionLevel _crvStepCompressionLevel;
  mu2e::CompressionLevel _simParticleCompressionLevel;
  mu2e::CompressionLevel _stepPointMCCompressionLevel;
  mu2e::CompressionLevel _genealogyCompressionLevel;
  mu2e::CompressionLevel _mcTrajectoryCompressionLevel;
  mu2e::CompressionLevel _primaryParticleCompressionLevel;

  // unique_ptrs to the new output collections
  std::unique_ptr<mu2e::StrawGasStepCollection> _newStrawGasSteps;
  std::unique_ptr<CaloShowerStepCollection> _newCaloShowerSteps;
  std::unique_ptr<CrvStepCollection> _newCrvSteps;
  std::map<InstanceLabel, std::unique_ptr<StepPointMCCollection> > _newStepPointMCs;
  std::unique_ptr<MCTrajectoryCollection> _newMCTrajectories;
  std::unique_ptr<PrimaryParticle> _newPrimaryParticle;

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
  SimParticleRemapping _simPtrRemap;
};


mu2e::CompressDetStepMCs::CompressDetStepMCs(const Parameters& conf)
  : art::EDProducer(conf),
    _conf(conf()),
    _strawGasStepCompressionLevel(mu2e::CompressionLevel::findByName(_conf.strawGasStepCompressionLevel())),
    _caloShowerStepCompressionLevel(mu2e::CompressionLevel::findByName(_conf.caloShowerStepCompressionLevel())),
    _crvStepCompressionLevel(mu2e::CompressionLevel::findByName(_conf.crvStepCompressionLevel())),
    _simParticleCompressionLevel(mu2e::CompressionLevel::findByName(_conf.simParticleCompressionLevel())),
  _stepPointMCCompressionLevel(mu2e::CompressionLevel::findByName(_conf.stepPointMCCompressionLevel())),
  _genealogyCompressionLevel(mu2e::CompressionLevel::findByName(_conf.genealogyCompressionLevel())),
  _mcTrajectoryCompressionLevel(mu2e::CompressionLevel::findByName(_conf.mcTrajectoryCompressionLevel())),
  _primaryParticleCompressionLevel(mu2e::CompressionLevel::findByName(_conf.primaryParticleCompressionLevel()))
{
  // Check that we have valid compression levels for this module
  checkCompressionLevels();

  // Call appropriate produces<>() functions here.
  produces<GenParticleCollection>();
  produces<SimParticleCollection>();

  produces<StrawGasStepCollection>();
  produces<CaloShowerStepCollection>();
  produces<CrvStepCollection>();

  for (const auto& i_tag : _conf.stepPointMCTags()) {
    produces<StepPointMCCollection>( i_tag.instance() );
  }
  produces<MCTrajectoryCollection>();
  produces<PrimaryParticle>();
}

void mu2e::CompressDetStepMCs::produce(art::Event & event)
{
  _newStrawGasSteps = std::unique_ptr<StrawGasStepCollection>(new StrawGasStepCollection);
  _newCaloShowerSteps = std::unique_ptr<CaloShowerStepCollection>(new CaloShowerStepCollection);
  _newCrvSteps = std::unique_ptr<CrvStepCollection>(new CrvStepCollection);

  for (const auto& i_tag : _conf.stepPointMCTags()) {
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
  _simPtrRemap.clear();

  // Compress detector steps and record which SimParticles we want to keep
  compressStrawGasSteps(event);
  compressCaloShowerSteps(event);
  compressCrvSteps(event);
  if (_stepPointMCCompressionLevel == mu2e::CompressionLevel::kNoCompression) {
    // if we are not compressing StepPointMCs, then
    // we want to make sure we record all their SimParticles
    compressStepPointMCs(event);
  }
  if (_mcTrajectoryCompressionLevel == mu2e::CompressionLevel::kNoCompression) {
    // if we are not compressing MCTrajectories, then
    // we want to make sure we record all their SimParticles
    compressMCTrajectories(event);
  }
  compressPrimaryParticle(event);

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
  if (_mcTrajectoryCompressionLevel == mu2e::CompressionLevel::kSimParticleCompression) {
    // if we are compressing MCTrajectories based on the SimParticles we are keeping,
    // then compress MCTrajectories now
    compressMCTrajectories(event);
  }

  // Update all the data products so that their SimParticlePtrs point to the new collection
  updateStrawGasSteps();
  updateCaloShowerSteps();
  updateCrvSteps();
  updateStepPointMCs();
  updateMCTrajectories();
  updatePrimaryParticle();

  // Now add everything to the event
  event.put(std::move(_newStrawGasSteps));
  event.put(std::move(_newCaloShowerSteps));
  event.put(std::move(_newCrvSteps));
  for (const auto& i_tag : _conf.stepPointMCTags()) {
    event.put(std::move(_newStepPointMCs.at(i_tag.instance())), i_tag.instance());
  }
  event.put(std::move(_newSimParticles));
  event.put(std::move(_newGenParticles));
  event.put(std::move(_newMCTrajectories));
  event.put(std::move(_newPrimaryParticle));
}

void mu2e::CompressDetStepMCs::compressStrawGasSteps(const art::Event& event) {
  art::Handle<mu2e::StrawGasStepCollection> strawGasStepsHandle;
  event.getByLabel(_conf.strawGasStepTag(), strawGasStepsHandle);
  const auto& strawGasSteps = *strawGasStepsHandle;
  if(_conf.debugLevel()>0 && strawGasSteps.size()>0) {
    std::cout << "Compressing StrawGasSteps from " << _conf.strawGasStepTag() << std::endl;
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
  art::Handle<mu2e::CaloShowerStepCollection> caloShowerStepsHandle;
  event.getByLabel(_conf.caloShowerStepTag(), caloShowerStepsHandle);
  const auto& caloShowerSteps = *caloShowerStepsHandle;
  if(_conf.debugLevel()>0 && caloShowerSteps.size()>0) {
    std::cout << "Compressing CaloShowerSteps from " << _conf.caloShowerStepTag() << std::endl;
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
  art::Handle<mu2e::CrvStepCollection> crvStepsHandle;
  event.getByLabel(_conf.crvStepTag(), crvStepsHandle);
  const auto& crvSteps = *crvStepsHandle;
  if(_conf.debugLevel()>0 && crvSteps.size()>0) {
    std::cout << "Compressing CrvSteps from " << _conf.crvStepTag() << std::endl;
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

void mu2e::CompressDetStepMCs::updateStrawGasSteps() {
  for (auto& i_strawGasStep : *_newStrawGasSteps) {
    const auto& oldSimPtr = i_strawGasStep.simParticle();
    art::Ptr<mu2e::SimParticle> newSimPtr = _simPtrRemap.at(oldSimPtr);
    if(_conf.debugLevel()>0) {
      std::cout << "Updating SimParticlePtr in StrawGasStep from " << oldSimPtr << " to " << newSimPtr << std::endl;
    }

    i_strawGasStep.simParticle() = newSimPtr;
  }
}

void mu2e::CompressDetStepMCs::updateCaloShowerSteps() {
  for (auto& i_caloShowerStep : *_newCaloShowerSteps) {
    const auto& oldSimPtr = i_caloShowerStep.simParticle();
    art::Ptr<mu2e::SimParticle> newSimPtr = _simPtrRemap.at(oldSimPtr);
    if(_conf.debugLevel()>0) {
      std::cout << "Updating SimParticlePtr in CaloShowerStep from " << oldSimPtr << " to " << newSimPtr << std::endl;
    }

    i_caloShowerStep.setSimParticle(newSimPtr);
  }
}

void mu2e::CompressDetStepMCs::updateCrvSteps() {
  for (auto& i_crvStep : *_newCrvSteps) {
    const auto& oldSimPtr = i_crvStep.simParticle();
    art::Ptr<mu2e::SimParticle> newSimPtr = _simPtrRemap.at(oldSimPtr);
    if(_conf.debugLevel()>0) {
      std::cout << "Updating SimParticlePtr in CrvStep from " << oldSimPtr << " to " << newSimPtr << std::endl;
    }

    i_crvStep.simParticle() = newSimPtr;
  }
}

void mu2e::CompressDetStepMCs::compressSimParticles(const art::Event& event) {
  // Now compress the SimParticleCollections into their new collections
  unsigned int keep_size = 0;
  const auto& oldSimParticles = event.getValidHandle<mu2e::SimParticleCollection>(_conf.simParticleTag());
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
  compressSimParticleCollection(_newSimParticlesPID, _newSimParticleGetter, *oldSimParticles, simPartSelector, *_newSimParticles);

  // Fill out the SimParticleRemapping
  for (const auto& i_keptSimPart : _simParticlesToKeep[i_product_id]) {
    cet::map_vector_key oldKey = cet::map_vector_key(i_keptSimPart.key());
    _simPtrRemap[i_keptSimPart] = art::Ptr<mu2e::SimParticle>(_newSimParticlesPID, oldKey.asUint(), _newSimParticleGetter);
    if (_conf.debugLevel()>0) {
      std::cout << "Compressing SimParticle " << i_keptSimPart << " --> " << _simPtrRemap[i_keptSimPart] << std::endl;
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

  // If we asked for the genealogy to be compressed, we will now end up with some missing links
  // add them back as truncated SimParticles
  if (_genealogyCompressionLevel != mu2e::CompressionLevel::kNoCompression) {
    // Go through the particles we are keeping and see if any parents are not there
    for (const auto& i_keptSimPart : _simParticlesToKeep[i_product_id]) {

      bool addedTruncated = false; // want to keep track of this because we don't want to add a truncated particle for every particle in the genealogy
      int n_truncated_added = 0;
      cet::map_vector_key truncKey;
      art::Ptr<mu2e::SimParticle> truncSimPtr;
      SimParticle* truncated = 0;

      art::Ptr<mu2e::SimParticle> i_childPtr = i_keptSimPart;//art::Ptr<mu2e::SimParticle>(_newSimParticlesPID, i_newSimParticle.first.asUint(), _newSimParticleGetter);
      art::Ptr<mu2e::SimParticle> i_parentPtr = i_childPtr->parent();
      while (i_parentPtr) {
        if (_simPtrRemap.find(i_parentPtr) == _simPtrRemap.end()) { // the parent is not in the output collection
          if (_conf.debugLevel()>0) {
            std::cout << "SimParticle " << i_parentPtr << " is not in output collection because it has been compressed away by genealogy compression" << std::endl;
          }

          // If we haven't added a truncated particle in this genealogy
          if (!addedTruncated) {
            // Add a truncated particle
            truncated = new SimParticle();
            truncKey = cet::map_vector_key(_conf.truncatedSimParticleKeyOffset() + n_truncated_added);
            truncated->id() = truncKey;
            //            truncated->addDaughter(i_childPtr);
            ++n_truncated_added;
            truncSimPtr = art::Ptr<mu2e::SimParticle>(_newSimParticlesPID, truncKey.asUint(), _newSimParticleGetter);
            if (_conf.debugLevel()>0) {
              std::cout << "Adding truncated SimParticle " << truncSimPtr << std::endl;
            }
            addedTruncated = true;
          }

          _simPtrRemap[i_parentPtr] = truncSimPtr;
          if (_conf.debugLevel()>0) {
            std::cout << "Previous SimParticlePtrs for " << i_parentPtr << " will now point to " << _simPtrRemap[i_parentPtr] << std::endl;
          }
        }
        else {
          // this parent is in the output collection so
          if (_conf.debugLevel()>0) {
            std::cout << "SimParticle " << i_parentPtr << " is in the output collection" << std::endl;
          }
          if(truncated) {
            truncated->parent() = i_parentPtr; // have the truncated particle's parent point to this
            addedTruncated = false; // and flag that we will need a new truncated particle next time
            (*_newSimParticles)[truncKey] = *truncated;

            if (_conf.debugLevel()>0) {
              std::cout << "Set truncated SimParticle (" << truncSimPtr << ") parent to " << truncated->parent() << std::endl;
            }
          }
        }
        i_parentPtr = i_parentPtr->parent();
      }
    }
  }

  // Fill out the SimParticleRemapping so that the removed links point to the truncated SimParticle
}

void mu2e::CompressDetStepMCs::compressGenParticles() {
  // Loop through the new SimParticles to keep any GenParticles
  for (auto& i_simParticle : *_newSimParticles) {
    mu2e::SimParticle& newsim = i_simParticle.second;
    if(newsim.genParticle().isNonnull()) { // will crash if not resolvable
      // Copy GenParticle to the new collection
      _newGenParticles->emplace_back(*newsim.genParticle());
      newsim.genParticle() = art::Ptr<mu2e::GenParticle>(_newGenParticlesPID, _newGenParticles->size()-1, _newGenParticleGetter);
    }
  }
}

void mu2e::CompressDetStepMCs::compressStepPointMCs(const art::Event& event) {

  for (const auto& i_tag : _conf.stepPointMCTags()) {
    const auto& stepPointMCs = event.getValidHandle<StepPointMCCollection>(i_tag);
    if(_conf.debugLevel()>0 && stepPointMCs->size()>0) {
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

  const auto& mcTrajectories = event.getValidHandle<MCTrajectoryCollection>(_conf.mcTrajectoryTag());
  if(_conf.debugLevel()>0 && mcTrajectories->size()>0) {
    std::cout << "Compressing MCTrajectories from " << _conf.mcTrajectoryTag() << std::endl;
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
            MCTrajectory newMCTrajectory(mcTrajectory.second);
            _newMCTrajectories->insert(std::pair<art::Ptr<SimParticle>, mu2e::MCTrajectory>(mcTrajectory.second.sim(), newMCTrajectory));
          }
        }
      }
    }
    else if (_mcTrajectoryCompressionLevel == mu2e::CompressionLevel::kNoCompression) {
      MCTrajectory newMCTrajectory(mcTrajectory.second);
      _newMCTrajectories->insert(std::pair<art::Ptr<SimParticle>, mu2e::MCTrajectory>(mcTrajectory.second.sim(), newMCTrajectory));
      recordSimParticle(mcTrajectory.second.sim());
    }
    else {
      throw cet::exception("CompressDetStepMCs") << "Unrecognized compression level \"" << _mcTrajectoryCompressionLevel.name() <<"\" for MCTrajectories" << std::endl;
    }
  }
}

void mu2e::CompressDetStepMCs::compressPrimaryParticle(const art::Event& event) {

  const auto& primaryParticle = event.getValidHandle<PrimaryParticle>(_conf.primaryParticleTag());
  if(_conf.debugLevel()>0) {
    std::cout << "Compressing PrimaryParticle from " << _conf.primaryParticleTag() << std::endl;
  }
  _newPrimaryParticle = std::unique_ptr<PrimaryParticle>(new PrimaryParticle(*primaryParticle));
  for (const auto& i_simParticle : _newPrimaryParticle->primarySimParticles()) {
    recordSimParticle(i_simParticle);
  }
}

void mu2e::CompressDetStepMCs::updateStepPointMCs() {
  for (const auto& i_tag : _conf.stepPointMCTags()) {
    for (auto& i_stepPointMC : *(_newStepPointMCs.at(i_tag.instance()))) {
      const auto& oldSimPtr = i_stepPointMC.simParticle();
      art::Ptr<mu2e::SimParticle> newSimPtr = _simPtrRemap.at(oldSimPtr);
      if(_conf.debugLevel()>0) {
        std::cout << "Updating SimParticlePtr in StepPointMC from " << oldSimPtr << " to " << newSimPtr << std::endl;
      }
      i_stepPointMC.simParticle() = newSimPtr;
    }
  }
}

void mu2e::CompressDetStepMCs::updateMCTrajectories() {
  for (auto& i_mcTrajectory : *_newMCTrajectories) {
    const auto& oldSimPtr = i_mcTrajectory.second.sim();
    art::Ptr<mu2e::SimParticle> newSimPtr = _simPtrRemap.at(oldSimPtr);
    if(_conf.debugLevel()>0) {
      std::cout << "Updating SimParticlePtr in MCTrajectory from " << oldSimPtr << " to " << newSimPtr << std::endl;
    }
    i_mcTrajectory.second.sim() = newSimPtr;
    auto mcTrajPair = _newMCTrajectories->extract(i_mcTrajectory.first);
    mcTrajPair.key() = newSimPtr;
    _newMCTrajectories->insert(std::move(mcTrajPair));
    //    i_mcTrajectory.first = newSimPtr;
  }
}

void mu2e::CompressDetStepMCs::updatePrimaryParticle() {
  for (auto& i_simParticlePtr : _newPrimaryParticle->modifySimParticles()) {
    art::Ptr<mu2e::SimParticle> newSimPtr = _simPtrRemap.at(i_simParticlePtr);
    if(_conf.debugLevel()>0) {
      std::cout << "Updating SimParticlePtr in PrimaryParticle from " << i_simParticlePtr << " to " << newSimPtr << std::endl;
    }
    i_simParticlePtr = newSimPtr;
  }
}

void mu2e::CompressDetStepMCs::recordSimParticle(const art::Ptr<mu2e::SimParticle>& sim_ptr) {
  // Also need to add all the parents too
  _simParticlesToKeep[sim_ptr.id()].insert(sim_ptr);
  art::Ptr<mu2e::SimParticle> childPtr = sim_ptr;
  art::Ptr<mu2e::SimParticle> parentPtr = childPtr->parent();

  if(_conf.debugLevel()>0) {
    std::cout << "Recording SimParticle " << sim_ptr << std::endl;
  }
  while (parentPtr) {
    if (_genealogyCompressionLevel == mu2e::CompressionLevel::kNoCompression) {
      _simParticlesToKeep[sim_ptr.id()].insert(parentPtr);
      if(_conf.debugLevel()>0) {
        std::cout << "and recording it's parent " << parentPtr << std::endl;
      }
    }
    else if (_genealogyCompressionLevel == mu2e::CompressionLevel::kFullCompression) {
      if (parentPtr->isPrimary()) {
        _simParticlesToKeep[sim_ptr.id()].insert(parentPtr);
        if(_conf.debugLevel()>0) {
          std::cout << "and recording it's parent " << parentPtr << std::endl;
        }
      }
      else {
        if(_conf.debugLevel()>0) {
          std::cout << "and *not* recording it's parent " << parentPtr << std::endl;
        }
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

  if (_simParticleCompressionLevel != mu2e::CompressionLevel::kNoCompression &&
      _simParticleCompressionLevel != mu2e::CompressionLevel::kFullCompression) {
    throw cet::exception("CompressDetStepMCs") << "This module does not allow SimParticles to be compressed with compression level \"" << _simParticleCompressionLevel.name() <<"\"" << std::endl;
  }

  if (_stepPointMCCompressionLevel != mu2e::CompressionLevel::kNoCompression &&
      _stepPointMCCompressionLevel != mu2e::CompressionLevel::kSimParticleCompression) {
    throw cet::exception("CompressDetStepMCs") << "This module does not allow StepPointMCs to be compressed with compression level \"" << _stepPointMCCompressionLevel.name() <<"\"" << std::endl;
  }

  if (_genealogyCompressionLevel != mu2e::CompressionLevel::kNoCompression &&
      _genealogyCompressionLevel != mu2e::CompressionLevel::kFullCompression) {
    throw cet::exception("CompressDetStepMCs") << "This module does not allow the genealogy to be compressed with compression level \"" << _genealogyCompressionLevel.name() <<"\"" << std::endl;
  }

  if (_mcTrajectoryCompressionLevel != mu2e::CompressionLevel::kNoCompression &&
      _mcTrajectoryCompressionLevel != mu2e::CompressionLevel::kSimParticleCompression) {
    throw cet::exception("CompressDetStepMCs") << "This module does not allow MCTrajectories to be compressed with compression level \"" << _mcTrajectoryCompressionLevel.name() <<"\"" << std::endl;
  }

  if (_primaryParticleCompressionLevel != mu2e::CompressionLevel::kNoCompression) {
    throw cet::exception("CompressDetStepMCs") << "This module does not allow PrimaryPartilces to be compressed with compression level \"" << _primaryParticleCompressionLevel.name() <<"\"" << std::endl;
  }

}

DEFINE_ART_MODULE(mu2e::CompressDetStepMCs)
