// Standalone ExtMonFNALGun generator.  Unlike the EventGenerator
// version it does not require a genconfig file and is instead
// configured via the framework, making parameter scans much easier to
// perform.
//
// Andrei Gaponenko, 2012

#include <string>
#include <memory>
#include <vector>

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "SeedService/inc/SeedService.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

#include "EventGenerator/inc/ExtMonFNALGunImpl.hh"

namespace mu2e {

  class ExtMonFNALGun : public art::EDProducer {
    std::vector<ExtMonFNALGunImpl::Config> conf_;

    // ExtMonFNALGun has a non-movable member ParticleGunImpl, thus hold them by pointers
    std::vector<std::unique_ptr<ExtMonFNALGunImpl> > guns_;
    CLHEP::HepRandomEngine& engine_;

    void produce(art::Event& event) override;
    void beginRun(art::Run& run) override;
  public:

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Sequence<fhicl::Table<ExtMonFNALGunImpl::Config> > guns {
        Name("guns"),
          Comment("A list of ExtMon guns to call.")
          };
    };

    using Parameters = art::EDProducer::Table<Config>;
    explicit ExtMonFNALGun(const Parameters& conf);
  };

  ExtMonFNALGun::ExtMonFNALGun(const Parameters& conf)
    : EDProducer{conf}
    , conf_(conf().guns())
    , engine_{createEngine(art::ServiceHandle<SeedService>{}->getSeed())}
    {
      produces<GenParticleCollection>();
      if(conf_.empty()) {
        throw cet::exception("BADCONFIG")<<"Error: no ExtMon guns defined.\n";
      }
    }

  void ExtMonFNALGun::beginRun(art::Run&) {
    guns_.reserve(conf_.size());
    for(const auto& c: conf_) {
      guns_.emplace_back(new mu2e::ExtMonFNALGunImpl{engine_, c});
    }
  }

  void ExtMonFNALGun::produce(art::Event& event) {
    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
    for(auto& gun : guns_) {
      gun->generate(*output);
    }
    event.put(std::move(output));
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNALGun);
