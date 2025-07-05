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
    // for applying selections to stops
    struct StopConfig {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::OptionalAtom<double> tmin{Name("tmin"), Comment("Minimum stop time")};
      fhicl::OptionalAtom<double> tmax{Name("tmax"), Comment("Maximum stop time")};
    };

    // top-level configuration
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag> simCollTag{Name("simParticles"),Comment("A SimParticleCollection with input stopped particles")};
      fhicl::Atom<int> pdgId{Name("pdgId"),Comment("Particle PDG code to select")};
      fhicl::Atom<int> processCode{Name("processCode"), Comment("Particle end process code to select")};
      fhicl::OptionalSequence<int> decayOffPdgs{Name("decayOffPdgs"), Comment("List of PDG IDs that decay was turned off, for event weights")};
      fhicl::OptionalTable<StopConfig> cuts { Name("cuts"), Comment("Stop selection table") };
      fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("Diagnostic print level"), 0};
    };

    // Cut selection struct
    struct StopCuts {
      using var=std::pair<bool, double>;
      var tmin_;
      var tmax_;
      StopCuts() : tmin_(var(false, 0.)), tmax_(var(false, 0.)) {}
      StopCuts(StopConfig config) : StopCuts(){
        tmin_.first = config.tmin(tmin_.second);
        tmax_.first = config.tmax(tmax_.second);
      }
      bool apply_cuts(art::Ptr<SimParticle> sim) {
        if(sim.isNull()) return false;
        const double t = sim->endGlobalTime();
        if(tmin_.first && t < tmin_.second) return false;
        if(tmax_.first && t > tmax_.second) return false;
        return true;
      }
    };

    explicit StopSelection(const art::EDProducer::Table<Config>& config);
    virtual void produce(art::Event& event) override;
    virtual void beginSubRun(art::SubRun& sr) override;
    virtual void endSubRun(art::SubRun& sr) override;
    virtual void beginJob() override;
    virtual void endJob() override;

  private:
    std::vector<art::Ptr<SimParticle>> pruneStops(const std::vector<art::Ptr<SimParticle>>& sims);

    art::ProductToken<SimParticleCollection> const simsToken_;
    PDGCode pdgId_;
    ProcessCode::enum_type processCode_;
    std::vector<int> decayOffPdgs_;
    int diagLevel_;
    StopCuts cuts_;
    SumOfWeights total_;
    int nempty_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat randomFlat_;
  };

  //--------------------------------------------------------------------------------------------------------------
  StopSelection::StopSelection(const art::EDProducer::Table<Config>& config) :
    EDProducer{config}
    , simsToken_{consumes<SimParticleCollection>(config().simCollTag())}
    , pdgId_{config().pdgId()}
    , processCode_{static_cast<ProcessCode::enum_type>(config().processCode())}
    , diagLevel_{config().diagLevel()}
    , nempty_(0)
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randomFlat_{eng_}
  {
    produces<mu2e::SimParticleCollection>();
    if(!config().decayOffPdgs(decayOffPdgs_)) decayOffPdgs_ = {};
    if(config().cuts()) cuts_ = StopCuts(*config().cuts());
    // if some decays were turned off, produce an event weight
    if(!decayOffPdgs_.empty()) produces<mu2e::EventWeight>();
    produces<SumOfWeights, art::InSubRun>();
  }

  //--------------------------------------------------------------------------------------------------------------
  void StopSelection::beginJob(){
  }

  //--------------------------------------------------------------------------------------------------------------
  // apply selections to the stop collection
  std::vector<art::Ptr<SimParticle>> StopSelection::pruneStops(const std::vector<art::Ptr<SimParticle>>& sims) {
    std::vector<art::Ptr<SimParticle>> outlist;
    for(auto sim : sims) {
      if(cuts_.apply_cuts(sim)) outlist.push_back(sim);
    }
    return outlist;
  }

  //--------------------------------------------------------------------------------------------------------------
  void StopSelection::produce(art::Event& event) {
    auto output{std::make_unique<SimParticleCollection>()};
    const bool make_weight = !decayOffPdgs_.empty();
    double weight = 1.;

    // get the sim collection and the corresponding stops
    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto stops = simParticleList(simh, pdgId_, processCode_);
    const auto selected = pruneStops(stops);

    if(diagLevel_ > 1) {
      printf("StopSelection::%s:: From %zu stops selected %zu candidates\n", __func__, stops.size(), selected.size());
    }
    // select a random stop if available
    if(!selected.empty()) {
      const int index = randomFlat_.fireInt(selected.size());
      art::Ptr<SimParticle> sim = selected[index];

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
      total_.add(weight);
    } else {
      ++nempty_;
    }


    // produce the data products
    event.put(std::move(output));
    if(make_weight) {
      std::unique_ptr<EventWeight> ew(new EventWeight(weight));
      event.put(std::move(ew));
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  void StopSelection::beginSubRun(art::SubRun&) {
    total_.reset();
  }

  //--------------------------------------------------------------------------------------------------------------
  void StopSelection::endSubRun(art::SubRun& sr) {
    sr.put(std::unique_ptr<SumOfWeights>(new SumOfWeights(total_.sum(), total_.count())), art::fullSubRun());
    if(diagLevel_ > 0) {
      printf("[StopSelection::%s] Selected %lu stops with a sum of weights of %.5g, %i empty events\n", __func__,
             total_.count(), total_.sum(), nempty_);
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  void StopSelection::endJob(){
  }
}

using mu2e::StopSelection;
DEFINE_ART_MODULE(StopSelection)
