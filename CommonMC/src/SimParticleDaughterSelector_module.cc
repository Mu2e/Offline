// Select daughers of a given set of particles based on their creation
// codes, and write the result into a new SimParticlePtrCollection.
//
// Andrei Gaponenko, 2014

#include <string>
#include <vector>
#include <utility>
#include <set>

#include "cetlib_except/exception.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

// Mu2e includes.
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/ProcessCode.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"

#include "TH1D.h"

namespace mu2e {
  namespace {
    std::string processToString(const SimParticle& p) {
      std::ostringstream os;

      if(static_cast<int>(p.pdgId()) < 1000000000) {
        os<<p.pdgId()<<" "<<ProcessCode(p.creationCode());
      }
      else {
        os<<"nucleus "<<ProcessCode(p.creationCode());
      }

      return os.str();
    }
  }

  class SimParticleDaughterSelector : public art::EDProducer {
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<art::InputTag> particleInput{ Name("particleInput"), Comment("The input collection.") };
      fhicl::Sequence<std::string> processes{ Name("processes"),
          Comment("A list of production process names, like \"DIO\", \"NuclearCapture\".\n"
                  "SimParticles with the given production codes will be written  to the output.\n"
                  "If the process list is empty, all particles are passed to the output.\n"
                  "You can use this mode to produce diagnostic histograms that will contain\nthe names of possible processes."
                  )
          };
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit SimParticleDaughterSelector(const Parameters& conf);

    void produce(art::Event& evt) override;
  private:
    art::InputTag particleInput_;
    std::set<ProcessCode::enum_type> procs_;

    art::ServiceHandle<art::TFileService> tfs() { return art::ServiceHandle<art::TFileService>(); }
    TH1* haccepted_;
    TH1* hignored_;
  };

  //================================================================
  SimParticleDaughterSelector::SimParticleDaughterSelector(const Parameters& conf)
    : EDProducer{conf}
    , particleInput_(conf().particleInput())
    , haccepted_(tfs()->make<TH1D>("accepted", "Accepted pdgId and process code pairs", 1, 0., 1.))
    , hignored_(tfs()->make<TH1D>("ignored", "Ignored pdgId and process code pairs", 1, 0., 1.))
  {
    produces<SimParticlePtrCollection>();
    const auto& processes = conf().processes();
    for(const auto& proc: processes) {
      procs_.insert(ProcessCode::findByName(proc).id());
    }

  }

  //================================================================
  void SimParticleDaughterSelector::produce(art::Event& event) {

    std::unique_ptr<SimParticlePtrCollection> output(new SimParticlePtrCollection());

    auto ih = event.getValidHandle<SimParticlePtrCollection>(particleInput_);
    for(const auto& i: *ih) {
      const auto& particle = *i;
      for(const auto& daughter: particle.daughters()) {
        const std::string strid(processToString(*daughter));
        if(procs_.empty() || (procs_.find(daughter->creationCode().id()) != procs_.end())) {
          output->emplace_back(daughter);
          haccepted_->Fill(strid.c_str(), 1.);
        }
        else {
          hignored_->Fill(strid.c_str(), 1.);
        }
      }
    }

    event.put(std::move(output));
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::SimParticleDaughterSelector);
