// Select a stopped particle and pass on a slimmed sim particle collection with this stop and its history
// author: M. MacKenzie  2025

// art includes
#include "art/Framework/Core/EDFilter.h"
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
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/SumOfWeights.hh"
#include "Offline/Mu2eUtilities/inc/SimParticleGetTau.hh"
#include "Offline/Mu2eUtilities/inc/simParticleList.hh"
#include "Offline/SeedService/inc/SeedService.hh"

// std includes
#include <iostream>
#include <string>

using namespace std;
namespace mu2e {

  class StopSelection : public art::EDFilter {
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
      fhicl::Atom<art::InputTag> stepCollTag{Name("stepPointMCs"),Comment("A StepPointMCCollection with input steps")};
      fhicl::Atom<int> pdgId{Name("pdgId"),Comment("Particle PDG code to select")};
      fhicl::Atom<int> processCode{Name("processCode"), Comment("Particle end process code to select")};
      fhicl::Atom<bool> filter{Name("filter"), Comment("Filter out events with no selected stops")};
      fhicl::OptionalSequence<int> decayOffPdgs{Name("decayOffPdgs"), Comment("List of PDG IDs that decay was turned off, for event weights")};
      fhicl::OptionalTable<StopConfig> cuts { Name("cuts"), Comment("Stop selection table") };
      fhicl::OptionalAtom<double> acceptRejectMax{Name("acceptRejectMax"), Comment("Perform accept/reject with a given maximum weight")};
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

    explicit StopSelection(const art::EDFilter::Table<Config>& config);
    virtual bool filter(art::Event& event) override;
    virtual bool beginSubRun(art::SubRun& sr) override;
    virtual bool endSubRun(art::SubRun& sr) override;
    virtual void beginJob() override;
    virtual void endJob() override;

  private:
    std::vector<art::Ptr<SimParticle>> pruneStops(const std::vector<art::Ptr<SimParticle>>& sims);
    void addSim(art::Ptr<SimParticle> sim, std::unique_ptr<SimParticleCollection>& coll,
                art::ProductID const& pid, art::EDProductGetter const* getter);
    void addSteps(const StepPointMCCollection& steps,
                  std::unique_ptr<SimParticleCollection>& sims,
                  std::unique_ptr<StepPointMCCollection>& out_steps,
                  art::ProductID const& pid, art::EDProductGetter const* getter);
    void cleanDaughters(std::unique_ptr<SimParticleCollection>& coll);

    art::ProductToken<SimParticleCollection> const simsToken_;
    art::ProductToken<StepPointMCCollection> const stepsToken_;
    PDGCode pdgId_;
    ProcessCode::enum_type processCode_;
    bool filter_;
    std::vector<int> decayOffPdgs_;
    double maxWeight_;
    bool acceptReject_;
    int diagLevel_;
    StopCuts cuts_;
    SumOfWeights total_;
    SumOfWeights sampled_;
    SumOfWeights selected_;
    int nempty_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    CLHEP::RandFlat randomFlat_;
  };

  //--------------------------------------------------------------------------------------------------------------
  StopSelection::StopSelection(const art::EDFilter::Table<Config>& config) :
    EDFilter{config}
    , simsToken_{consumes<SimParticleCollection>(config().simCollTag())}
    , stepsToken_{consumes<StepPointMCCollection>(config().stepCollTag())}
    , pdgId_{config().pdgId()}
    , processCode_{static_cast<ProcessCode::enum_type>(config().processCode())}
    , filter_{config().filter()}
    , acceptReject_{config().acceptRejectMax(maxWeight_)}
    , diagLevel_{config().diagLevel()}
    , nempty_(0)
    , eng_{createEngine(art::ServiceHandle<SeedService>()->getSeed())}
    , randomFlat_{eng_}
  {
    produces<mu2e::SimParticleCollection>();
    produces<mu2e::StepPointMCCollection>();
    if(!config().decayOffPdgs(decayOffPdgs_)) decayOffPdgs_ = {};
    if(config().cuts()) cuts_ = StopCuts(*config().cuts());
    // if some decays were turned off, produce an event weight
    if(!decayOffPdgs_.empty()) produces<mu2e::EventWeight>();
    produces<SumOfWeights, art::InSubRun>("total");
    produces<SumOfWeights, art::InSubRun>("sampled");
    if(diagLevel_ > 0 && acceptReject_) {
      printf("[StopSelection::%s] Performing accept/rejection selection with a maximum weight of %.3g\n", __func__, maxWeight_);
    }
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
  // add a particle and its lineage to the collection
  void StopSelection::addSim(art::Ptr<SimParticle> sim, std::unique_ptr<SimParticleCollection>& coll,
                             art::ProductID const& pid, art::EDProductGetter const* getter) {
      while(sim.isNonnull()) {
        // check if the sim has already been added
        if(coll->has(sim->id())) break;
        //if not, add the sim to the output and continue up the lineage
        auto sim_id_pair = std::make_pair(sim->id(), *sim);
        auto& new_sim = sim_id_pair.second;
        auto parent = sim->parent();
        new_sim.parent() = art::Ptr<SimParticle>(pid, parent.key(), getter);
        for (auto& daughter : new_sim.daughters()) {
          daughter = art::Ptr<SimParticle>(pid,daughter.key(),getter);
        }
        sim = parent;
        coll->insert(sim_id_pair);
      }
  }


  //--------------------------------------------------------------------------------------------------------------
  // add relevant step points
  void StopSelection::addSteps(const StepPointMCCollection& steps,
                               std::unique_ptr<SimParticleCollection>& sims,
                               std::unique_ptr<StepPointMCCollection>& out_steps,
                               art::ProductID const& pid, art::EDProductGetter const* getter) {
    // for every input step point, add it to the output if relevant, with remapping
    for(auto& step : steps) {
      bool found = false;
      const auto id = step.simParticle()->id();
      for(auto& sim : *sims) {
        if(sim.first == id) {
          found = true;
          break;
        }
      }
      if(found) { // save the step
        auto out_step(step);
        out_step.simParticle() = art::Ptr<SimParticle>(pid, step.simParticle().key(), getter);
        out_steps->push_back(out_step);
      }
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  // remove daughters that were not added to the collection
  void StopSelection::cleanDaughters(std::unique_ptr<SimParticleCollection>& coll) {
    for(auto& sim : *coll) {
      std::vector<art::Ptr<SimParticle>> out_daughters;
      for (auto& daughter : sim.second.daughters()) {
        if(coll->has(cet::map_vector_key(daughter.key()))) out_daughters.push_back(daughter);
      }
      sim.second.setDaughterPtrs(out_daughters);
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  bool StopSelection::filter(art::Event& event) {
    auto out_sims{std::make_unique<SimParticleCollection>()};
    auto out_steps{std::make_unique<StepPointMCCollection>()};
    const PhysicsParams& gc = *GlobalConstantsHandle<PhysicsParams>();
    const bool make_weight = !decayOffPdgs_.empty();
    double weight = 1.;

    // for mapping new sim art Ptrs
    auto out_pid = event.getProductID<SimParticleCollection>("");
    auto out_getter = event.productGetter(out_pid);

    // get the sim collection and the corresponding stops
    const auto simh = event.getValidHandle<SimParticleCollection>(simsToken_);
    const auto steps = event.getValidHandle<StepPointMCCollection>(stepsToken_);
    const auto stops = simParticleList(simh, pdgId_, processCode_); // stopped particles
    const auto selected = pruneStops(stops); // stopped particles satisfying a given selection

    if(diagLevel_ > 1) {
      printf("[StopSelection::%s] From %zu stops selected %zu candidates\n", __func__, stops.size(), selected.size());
    }

    const size_t nstops = selected.size();
    const int offset = (acceptReject_) ? 0 : randomFlat_.fireInt(nstops); // if picking 1 stop, make it random, otherwise sample all
    int naccepted = 0;
    for(size_t istop = 0; istop < nstops; ++istop) {
      const size_t index = (istop + offset) % nstops;
      art::Ptr<SimParticle> sim = selected[index];

      // get the selected stop's time weight if requested
      if(make_weight) {
        const double tau = SimParticleGetTau::calculate(sim, decayOffPdgs_, gc);
        weight = std::exp(-tau);
      }
      total_.add(weight);
      if(istop > 0 && !acceptReject_) continue; // if not doing accept reject, don't need to test more, just store their weights
      sampled_.add(weight);

      // check if doing accept/reject
      if(acceptReject_) {
        if(weight < maxWeight_ * randomFlat_.fire()) continue; // if it fails, continue to the next stop
        weight = 1.; // reset as no event weight is used if using accept/reject
      }
      ++naccepted;

      // save the entire lineage, updating the pointer mapping
      addSim(sim, out_sims, out_pid, out_getter);

      selected_.add(weight);
    }
    if(nstops == 0) ++nempty_;

    // for every input step point, add it to the output if relevant, with remapping
    if(out_sims->size() > 0) {
      addSteps(*steps, out_sims, out_steps, out_pid, out_getter);
      cleanDaughters(out_sims);
    }

    if(diagLevel_ > 1) {
      printf("[StopSelection::%s] From %zu candidates accepted %i stops\n", __func__, selected.size(), naccepted);
    }

    // determine whether or not to accept the event
    const bool pass = !filter_ || out_sims->size() > 0;

    // produce the data products
    event.put(std::move(out_sims));
    event.put(std::move(out_steps));
    if(make_weight) {
      event.put(std::unique_ptr<EventWeight>(new EventWeight(weight)));
    }

    return pass;
  }

  //--------------------------------------------------------------------------------------------------------------
  bool StopSelection::beginSubRun(art::SubRun&) {
    total_.reset();
    selected_.reset();
    return true;
  }

  //--------------------------------------------------------------------------------------------------------------
  bool StopSelection::endSubRun(art::SubRun& sr) {
    sr.put(std::unique_ptr<SumOfWeights>(new SumOfWeights(total_.sum(), total_.count())), "total", art::fullSubRun());
    sr.put(std::unique_ptr<SumOfWeights>(new SumOfWeights(sampled_.sum(), sampled_.count())), "sampled", art::fullSubRun());
    if(diagLevel_ > 0) {
      printf("[StopSelection::%s] Selected %lu stops with a sum of weights of %.5g, %lu total events with %.5g sum of weights, sampled %lu events with %.5g sum of weights, %i empty events\n", __func__,
             selected_.count(), selected_.sum(), total_.count(), total_.sum(), sampled_.count(), sampled_.sum(), nempty_);
    }
    return true;
  }

  //--------------------------------------------------------------------------------------------------------------
  void StopSelection::endJob(){
  }
}

using mu2e::StopSelection;
DEFINE_ART_MODULE(StopSelection)
