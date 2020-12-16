//
// Compression module that takes detector steps (e.g. StrawGasStep)
// and compresses out unwanted MC information (e.g. SimParticles)
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
}


class mu2e::CompressDetStepMCs : public art::EDProducer {
public:
  struct Config {
    using Name=fhicl::Name;
    using Comment=fhicl::Comment;

    // fhicl parameters for all of our inputs
    fhicl::Atom<art::InputTag> strawGasStepTag{Name("strawGasStepTag"), Comment("InputTag for the StrawGasSteps")};
    //    fhicl::Atom<art::InputTag> caloShowerStepTag{Name("caloShowerStepTag"), Comment("InputTag for CaloShowerSteps")};
    //    fhicl::Atom<art::InputTag> crvBarStepTag{Name("crvBarStepTag"), Comment("InputTag for CrvBarSteps")};
    fhicl::Atom<art::InputTag> simParticleTag{Name("simParticleTag"), Comment("InputTag for the SimParticleCollection")};
    fhicl::Atom<int> debugLevel{Name("debugLevel"), Comment("Debug level (0 = no debug output)")};
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
  void compressSimParticles(const art::Event& event);
  void compressGenParticles();
  void recordSimParticle(const art::Ptr<mu2e::SimParticle>& sim_ptr);

private:

  Config _conf;

  // unique_ptrs to the new output collections
  std::unique_ptr<mu2e::StrawGasStepCollection> _newStrawGasSteps;
  //  std::unique_ptr<CaloShowerStepCollection> _newCaloShowerSteps;
  //  std::unique_ptr<CrvBarStepCollection> _newCrvBarSteps;

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
    _conf(conf())
{
  // Call appropriate produces<>() functions here.
  produces<GenParticleCollection>();
  produces<SimParticleCollection>();

  produces<StrawGasStepCollection>();
  //  produces<CaloShowerStepCollection>();
  //  produces<CrvBarStepCollection>();
}

void mu2e::CompressDetStepMCs::produce(art::Event & event)
{
  _newStrawGasSteps = std::unique_ptr<StrawGasStepCollection>(new StrawGasStepCollection);
  //  _newCaloShowerSteps = std::unique_ptr<CaloShowerStepCollection>(new CaloShowerStepCollection);
  //  _newCrvBarSteps = std::unique_ptr<CrvBarStepCollection>(new CrvBarStepCollection);
  //  _newStrawGasStepsPID = event.getProductID<StrawGasStepCollection>();
//  _newStrawGasStepGetter = event.productGetter(_newStrawGasStepsPID);

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

  // Compress the SimParticles and record their new keys
  compressSimParticles(event);

  // Create the new GenParticleCollection for the SimParticles we are keeping
  // and have the new SimParticle point to the new GenParticle
  compressGenParticles();

  // Update all the detector steps so that their SimParticlePtrs point to the new collection
  updateStrawGasSteps();

  // Now add everything to the event
  event.put(std::move(_newStrawGasSteps));
  event.put(std::move(_newSimParticles));
  event.put(std::move(_newGenParticles));
}

void mu2e::CompressDetStepMCs::compressStrawGasSteps(const art::Event& event) {
  art::Handle<mu2e::StrawGasStepCollection> strawGasStepsHandle;
  event.getByLabel(_conf.strawGasStepTag(), strawGasStepsHandle);
  const auto& strawGasSteps = *strawGasStepsHandle;
  if(_conf.debugLevel()>0 && strawGasSteps.size()>0) {
    std::cout << "Compressing StrawGasSteps from " << _conf.strawGasStepTag() << std::endl;
  }
  for (const auto& i_strawGasStep : strawGasSteps) {
    recordSimParticle(i_strawGasStep.simParticle());
    StrawGasStep newStrawGasStep(i_strawGasStep);
    _newStrawGasSteps->push_back(newStrawGasStep);
  }
  if (_newStrawGasSteps->size() != strawGasSteps.size()) {
    throw cet::exception("CompressDetStepMCs") << "Number of StrawGasSteps in output collection (" << _newStrawGasSteps->size() << ") does not match the number of StrawGasSteps in the inout collection (" << strawGasSteps.size() << ")" << std::endl;
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

void mu2e::CompressDetStepMCs::compressSimParticles(const art::Event& event) {
  // Now compress the SimParticleCollections into their new collections
  unsigned int keep_size = 0;
  const auto& oldSimParticles = event.getValidHandle<mu2e::SimParticleCollection>(_conf.simParticleTag());
  art::ProductID i_product_id = oldSimParticles.id();
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


DEFINE_ART_MODULE(mu2e::CompressDetStepMCs)
