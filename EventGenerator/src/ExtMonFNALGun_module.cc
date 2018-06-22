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

#include "EventGenerator/inc/ExtMonFNALGun.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    class ExtMonFNALGun : public art::EDProducer {
      fhicl::ParameterSet pset_;

      // ExtMonFNALGun has a non-movable member ParticleGunImpl, thus hold them by pointers
      typedef std::vector<std::unique_ptr<mu2e::ExtMonFNALGun> > PGuns;
      PGuns guns_;
      CLHEP::HepRandomEngine& engine_;

      void produce(art::Event& event) override;
      void beginRun(art::Run& run) override;
    public:
      explicit ExtMonFNALGun(fhicl::ParameterSet const& pset);
    };

    ExtMonFNALGun::ExtMonFNALGun(fhicl::ParameterSet const& pset)
      : pset_(pset)
      , engine_{createEngine(art::ServiceHandle<SeedService>{}->getSeed())}
    {
      produces<GenParticleCollection>();
    }

    void ExtMonFNALGun::beginRun(art::Run&) {
      typedef std::vector<fhicl::ParameterSet> VGPars;
      VGPars vgp(pset_.get<VGPars>("guns"));

      guns_.reserve(vgp.size());
      for(unsigned i=0; i<vgp.size(); ++i) {
        guns_.emplace_back(new mu2e::ExtMonFNALGun{engine_, vgp[i]});
      }
    }

    void ExtMonFNALGun::produce(art::Event& event) {
      std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);
      for(auto& gun : guns_) {
        gun->generate(*output);
      }
      event.put(std::move(output));
    }

  } // namespace ExtMonFNAL
} // namespace mu2e

DEFINE_ART_MODULE(mu2e::ExtMonFNAL::ExtMonFNALGun);
