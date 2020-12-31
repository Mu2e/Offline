//
// Compression module that takes detector steps (StrawGasSteps, CaloShowerSteps, and CrvBarSteps)
// and compresses out unwanted MC information (StepPointMCs, SimParticles and GenParticles).
// Note that the number of detector steps is NOT reduced.
//
// Allowed compression levels for different data products:
// - StrawGasSteps : noCompression
// - CaloShowerSteps : noCompression
// - CrvBarSteps : noCompression
// - SimParticles : noCompression, fullCompression
// - StepPointMCs : noCompression, simParticleCompression
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
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "Mu2eUtilities/inc/compressSimParticleCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "Compression/inc/CompressionLevel.hh"

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

  typedef mu2e::StepPointMC CrvBarStep; // this will become a CrvBarStep class
  typedef mu2e::StepPointMCCollection CrvBarStepCollection;
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
    fhicl::Atom<art::InputTag> crvBarStepTag{Name("crvBarStepTag"), Comment("InputTag for CrvBarSteps")};
    fhicl::Sequence<art::InputTag> stepPointMCTags{Name("stepPointMCTags"), Comment("Sequence of InputTags for StepPointMCCollections (e.g. virtualdetector)")};
    fhicl::Atom<art::InputTag> simParticleTag{Name("simParticleTag"), Comment("InputTag for the SimParticleCollection")};
    fhicl::Atom<int> debugLevel{Name("debugLevel"), Comment("Debug level (0 = no debug output)")};
    fhicl::Atom<std::string> strawGasStepCompressionLevel{Name("strawGasStepCompressionLevel"), Comment("Compression level for StrawGasSteps")};
    fhicl::Atom<std::string> caloShowerStepCompressionLevel{Name("caloShowerStepCompressionLevel"), Comment("Compression level for CaloShowerSteps")};
    fhicl::Atom<std::string> crvBarStepCompressionLevel{Name("crvBarStepCompressionLevel"), Comment("Compression level for CrvBarSteps")};
    fhicl::Atom<std::string> simParticleCompressionLevel{Name("simParticleCompressionLevel"), Comment("Compression level for SimParticles")};
    fhicl::Atom<std::string> stepPointMCCompressionLevel{Name("stepPointMCCompressionLevel"), Comment("Compression level for StepPointMCs")};
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
  void compressCrvBarSteps(const art::Event& event);
  void updateCrvBarSteps();
  void compressStepPointMCs(const art::Event& event);
  void updateStepPointMCs();
  void compressSimParticles(const art::Event& event);
  void compressGenParticles();
  void recordSimParticle(const art::Ptr<mu2e::SimParticle>& sim_ptr);
  void checkCompressionLevels();

private:

  Config _conf;
  mu2e::CompressionLevel _strawGasStepCompressionLevel;
  mu2e::CompressionLevel _caloShowerStepCompressionLevel;
  mu2e::CompressionLevel _crvBarStepCompressionLevel;
  mu2e::CompressionLevel _simParticleCompressionLevel;
  mu2e::CompressionLevel _stepPointMCCompressionLevel;

  // unique_ptrs to the new output collections
  std::unique_ptr<mu2e::StrawGasStepCollection> _newStrawGasSteps;
  std::unique_ptr<CaloShowerStepCollection> _newCaloShowerSteps;
  std::unique_ptr<CrvBarStepCollection> _newCrvBarSteps;
  std::map<InstanceLabel, std::unique_ptr<StepPointMCCollection> > _newStepPointMCs;

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
    _crvBarStepCompressionLevel(mu2e::CompressionLevel::findByName(_conf.crvBarStepCompressionLevel())),
    _simParticleCompressionLevel(mu2e::CompressionLevel::findByName(_conf.simParticleCompressionLevel())),
    _stepPointMCCompressionLevel(mu2e::CompressionLevel::findByName(_conf.stepPointMCCompressionLevel()))
{
  // Check that we have valid compression levels for this module
  checkCompressionLevels();

  // Call appropriate produces<>() functions here.
  produces<GenParticleCollection>();
  produces<SimParticleCollection>();

  produces<StrawGasStepCollection>();
  produces<CaloShowerStepCollection>();
  produces<CrvBarStepCollection>("CRV"); // need to give an instance name because this is currently a StepPointMC

  for (const auto& i_tag : _conf.stepPointMCTags()) {
    produces<StepPointMCCollection>( i_tag.instance() );
  }
}

void mu2e::CompressDetStepMCs::produce(art::Event & event)
{
  _newStrawGasSteps = std::unique_ptr<StrawGasStepCollection>(new StrawGasStepCollection);
  _newCaloShowerSteps = std::unique_ptr<CaloShowerStepCollection>(new CaloShowerStepCollection);
  _newCrvBarSteps = std::unique_ptr<CrvBarStepCollection>(new CrvBarStepCollection);

  for (const auto& i_tag : _conf.stepPointMCTags()) {
    _newStepPointMCs[i_tag.instance()] = std::unique_ptr<StepPointMCCollection>(new StepPointMCCollection);
  }

  _newSimParticles = std::unique_ptr<SimParticleCollection>(new SimParticleCollection);
  _newSimParticlesPID = event.getProductID<SimParticleCollection>();
  _newSimParticleGetter = event.productGetter(_newSimParticlesPID);

  _newGenParticles = std::unique_ptr<GenParticleCollection>(new GenParticleCollection);
  _newGenParticlesPID = event.getProductID<GenParticleCollection>();
  _newGenParticleGetter = event.productGetter(_newGenParticlesPID);

  _simParticlesToKeep.clear();
  _simPtrRemap.clear();

  // Compress detector steps and record which SimParticles we want to keep
  compressStrawGasSteps(event);
  compressCaloShowerSteps(event);
  compressCrvBarSteps(event);
  if (_stepPointMCCompressionLevel == mu2e::CompressionLevel::kNoCompression) {
    // if we are not compressing StepPointMCs, then
    // we want to make sure we record all their SimParticles
    compressStepPointMCs(event);
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

  // Update all the data products so that their SimParticlePtrs point to the new collection
  updateStrawGasSteps();
  updateCaloShowerSteps();
  updateCrvBarSteps();
  updateStepPointMCs();

  // Now add everything to the event
  event.put(std::move(_newStrawGasSteps));
  event.put(std::move(_newCaloShowerSteps));
  event.put(std::move(_newCrvBarSteps), "CRV");
  for (const auto& i_tag : _conf.stepPointMCTags()) {
    event.put(std::move(_newStepPointMCs.at(i_tag.instance())), i_tag.instance());
  }
  event.put(std::move(_newSimParticles));
  event.put(std::move(_newGenParticles));
}

void mu2e::CompressDetStepMCs::compressStrawGasSteps(const art::Event& event) {
  art::Handle<mu2e::StrawGasStepCollection> strawGasStepsHandle;
  event.getByLabel(_conf.strawGasStepTag(), strawGasStepsHandle);
  if (strawGasStepsHandle.isValid()) {
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
}

void mu2e::CompressDetStepMCs::compressCaloShowerSteps(const art::Event& event) {
  art::Handle<mu2e::CaloShowerStepCollection> caloShowerStepsHandle;
  event.getByLabel(_conf.caloShowerStepTag(), caloShowerStepsHandle);
  if (caloShowerStepsHandle.isValid()) {
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
}

void mu2e::CompressDetStepMCs::compressCrvBarSteps(const art::Event& event) {
  art::Handle<mu2e::CrvBarStepCollection> crvBarStepsHandle;
  event.getByLabel(_conf.crvBarStepTag(), crvBarStepsHandle);
  const auto& crvBarSteps = *crvBarStepsHandle;
  if (crvBarStepsHandle.isValid()) {
    if(_conf.debugLevel()>0 && crvBarSteps.size()>0) {
      std::cout << "Compressing CrvBarSteps from " << _conf.crvBarStepTag() << std::endl;
    }
    for (const auto& i_crvBarStep : crvBarSteps) {
      if (_simParticleCompressionLevel == mu2e::CompressionLevel::kFullCompression) {
        recordSimParticle(i_crvBarStep.simParticle());
      }
      CrvBarStep newCrvBarStep(i_crvBarStep);
      _newCrvBarSteps->push_back(newCrvBarStep);
    }
    if(_crvBarStepCompressionLevel == mu2e::CompressionLevel::kNoCompression) {
      if (_newCrvBarSteps->size() != crvBarSteps.size()) {
        throw cet::exception("CompressDetStepMCs") << "Number of CrvBarSteps in output collection (" << _newCrvBarSteps->size() << ") does not match the number of CrvBarSteps in the input collection (" << crvBarSteps.size() << ") even though no compression has been requested (crvBarStepCompressionLevel = \"" << _crvBarStepCompressionLevel.name() << "\")" << std::endl;
      }
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

void mu2e::CompressDetStepMCs::updateCrvBarSteps() {
  for (auto& i_crvBarStep : *_newCrvBarSteps) {
    const auto& oldSimPtr = i_crvBarStep.simParticle();
    art::Ptr<mu2e::SimParticle> newSimPtr = _simPtrRemap.at(oldSimPtr);
    if(_conf.debugLevel()>0) {
      std::cout << "Updating SimParticlePtr in CrvBarStep from " << oldSimPtr << " to " << newSimPtr << std::endl;
    }

    i_crvBarStep.simParticle() = newSimPtr;
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

void mu2e::CompressDetStepMCs::recordSimParticle(const art::Ptr<mu2e::SimParticle>& sim_ptr) {
  // Also need to add all the parents too
  _simParticlesToKeep[sim_ptr.id()].insert(sim_ptr);
  art::Ptr<mu2e::SimParticle> childPtr = sim_ptr;
  art::Ptr<mu2e::SimParticle> parentPtr = childPtr->parent();

  if(_conf.debugLevel()>0) {
    std::cout << "Recording SimParticle " << sim_ptr << std::endl;
  }
  while (parentPtr) {
    if(_conf.debugLevel()>0) {
      std::cout << "and it's parent " << parentPtr << std::endl;
    }
    _simParticlesToKeep[sim_ptr.id()].insert(parentPtr);
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

  if (_crvBarStepCompressionLevel != mu2e::CompressionLevel::kNoCompression) {
    throw cet::exception("CompressDetStepMCs") << "This module does not allow CrvBarSteps to be compressed with compression level \"" << _crvBarStepCompressionLevel.name() <<"\"" << std::endl;
  }

  if (_simParticleCompressionLevel != mu2e::CompressionLevel::kNoCompression &&
      _simParticleCompressionLevel != mu2e::CompressionLevel::kFullCompression) {
    throw cet::exception("CompressDetStepMCs") << "This module does not allow SimParticles to be compressed with compression level \"" << _simParticleCompressionLevel.name() <<"\"" << std::endl;
  }

  if (_stepPointMCCompressionLevel != mu2e::CompressionLevel::kNoCompression &&
      _stepPointMCCompressionLevel != mu2e::CompressionLevel::kSimParticleCompression) {
    throw cet::exception("CompressDetStepMCs") << "This module does not allow StepPointMCs to be compressed with compression level \"" << _stepPointMCCompressionLevel.name() <<"\"" << std::endl;
  }
}

DEFINE_ART_MODULE(mu2e::CompressDetStepMCs)
