// Select a stopped particle and pass on a slimmed sim particle collection with this stop and its history
// author: M. MacKenzie  2025

// art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "art_root_io/TFileService.h"

// CLHEP
#include "CLHEP/Random/RandFlat.h"

// Mu2e includes
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/MCDataProducts/inc/EventWeight.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/SumOfWeights.hh"
#include "Offline/Mu2eUtilities/inc/SimParticleGetTau.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"
#include "Offline/SeedService/inc/SeedService.hh"

// std includes
#include <iostream>
#include <string>

using namespace std;
namespace mu2e {

  class StopSelection : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag> simCollTag{Name("simParticles"),Comment("A SimParticleCollection with input stopped particles")};
      fhicl::Atom<int> pdgId{Name("pdgId"),Comment("Particle PDG code to select")};
      fhicl::Atom<int> processCode{Name("processCode"), Comment("Particle end process code to select")};
      fhicl::OptionalSequence<int> decayOffPdgs{Name("decayOffPdgs"), Comment("List of PDG IDs that decay was turned off, for event weights")};
      fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("Diagnostic print level"), 0};
    };
    explicit StopSelection(const art::EDProducer::Table<Config>& config);
    virtual void produce(art::Event& event) override;
    virtual void beginSubRun(art::SubRun& sr) override;
    virtual void endSubRun(art::SubRun& sr) override;
    virtual void beginJob() override;
    virtual void endJob() override;

  private:
    art::ProductToken<SimParticleCollection> const simsToken_;
    PDGCode pdgId_;
    ProcessCode::enum_type processCode_;
    std::vector<int> decayOffPdgs_;
    int diagLevel_;
    SumOfWeights total_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat randomFlat_;
  };

  StopSelection::StopSelection(const art::EDProducer::Table<Config>& config) :
    EDProducer{config}
    , simsToken_{consumes<SimParticleCollection>(config().simCollTag())}
    , pdgId_{config().pdgId()}
    , processCode_{static_cast<ProcessCode::enum_type>(config().processCode())}
    , diagLevel_{config().diagLevel()}
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randomFlat_{eng_}
  {
    produces<mu2e::SimParticleCollection>();
    if(!config().decayOffPdgs(decayOffPdgs_)) decayOffPdgs_ = {};
    // if some decays were turned off, produce an event weight
    if(!decayOffPdgs_.empty()) produces<mu2e::EventWeight>();
    produces<SumOfWeights, art::InSubRun>();
  }

  void StopSelection::beginJob(){
  }

  void StopSelection::produce(art::Event& event) {
    auto output{std::make_unique<SimParticleCollection>()};
    const bool make_weight = !decayOffPdgs_.empty();
    double weight = 1.;

    // get the sim collection and the corresponding stops
    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto stops = simParticleList(simh, pdgId_, processCode_);

    if(diagLevel_ > 1) {
      printf("StopSelection::%s:: Printing input sim collections:\n", __func__);
    }
    // select a random stop if available
    if(!stops.empty()) {
      const int index = randomFlat_.fireInt(stops.size());
      art::Ptr<SimParticle> sim = stops[index];

      // get the selected stop's time weight if requested
      if(make_weight) {
        const PhysicsParams& gc = *GlobalConstantsHandle<PhysicsParams>();
        const double tau = SimParticleGetTau::calculate(sim, decayOffPdgs_, gc);
        weight = std::exp(-tau);
      }

      // save the entire lineage
      while(sim.isNonnull()) {
        // output->insert(std::make_pair(output->size(), *sim));
        output->insert(std::make_pair(sim->id(), *sim));
        sim = sim->parent();
      }
    }

    total_.add(weight);

    // produce the data products
    event.put(std::move(output));
    if(make_weight) {
      std::unique_ptr<EventWeight> ew(new EventWeight(weight));
      event.put(std::move(ew));
    }
  }

  void StopSelection::beginSubRun(art::SubRun&) {
    total_.reset();
  }

  void StopSelection::endSubRun(art::SubRun& sr) {
    sr.put(std::unique_ptr<SumOfWeights>(new SumOfWeights(total_.sum(), total_.count())), art::fullSubRun());
    if(diagLevel_ > 0) {
      printf("[StopSelection::%s] Selected %lu stops with a sum of weights of %.5g\n", __func__, total_.count(), total_.sum());
    }
  }

  void StopSelection::endJob(){
  }
}

using mu2e::StopSelection;
DEFINE_ART_MODULE(StopSelection)
