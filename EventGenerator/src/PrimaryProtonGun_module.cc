
/*
  This was a Replicated Module and was reverted to Legacy on Mar 27, 2022
  until GitHub issues #744 and #745 can be resolved.

  A plug-in for running PrimaryProtonGun-based event generator for running in MT art.
  It produces a GenParticleCollection of primary protons using the PrimaryProtonGun.

  These Collections are used in Mu2eG4_module.cc.

 Original author: Lisa Goodenough
  Date: 2020/02/19

*/

// Mu2e includes.
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/MCDataProducts/inc/GenParticle.hh"

// Particular generators that this code knows about.
#include "Offline/EventGenerator/inc/PrimaryProtonGunImpl.hh"
#include "Offline/SeedService/inc/SeedService.hh"

// Includes from art and its toolchain.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes.
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

namespace mu2e {

  class PrimaryProtonGun : public art::EDProducer {

  public:

      typedef art::EDProducer::Table<PrimaryProtonGunImpl::PrimaryProtonGunConfig> Parameters;

      explicit PrimaryProtonGun(Parameters const& params );
      // Accept compiler written d'tor.  Modules are never moved or copied.

      virtual void produce(art::Event& e ) override;
      virtual void beginRun(art::Run& r ) override;

  private:

      const PrimaryProtonGunImpl::PrimaryProtonGunConfig _config;

      CLHEP::HepRandomEngine& _engine;
      std::unique_ptr<PrimaryProtonGunImpl> _primaryProtonGunGenerator;

      // Number of times BeginRun is called on this module
      int ncalls = 0;

  };

    PrimaryProtonGun::PrimaryProtonGun(Parameters const& params ):
        art::EDProducer{params},
        _config(params()),
        _engine{createEngine(art::ServiceHandle<SeedService>{}->getSeed())}
    {
        produces<GenParticleCollection>();
    }

    void PrimaryProtonGun::beginRun(art::Run& run ){

    // The configuration of the PPG Generator does not change within a job.
    if ( ++ncalls > 1){
      mf::LogInfo("PrimaryProtonGun")
        << ", PrimaryProtonGun Generator does not change state at beginRun.  Hope that's OK.";
      return;
    }

    // Instantiate generator for this run.
    _primaryProtonGunGenerator = std::make_unique <PrimaryProtonGunImpl>( _engine, _config);

  }//beginRun

    void PrimaryProtonGun::produce(art::Event& evt ) {

        // Make the collections to hold the output.
        unique_ptr<GenParticleCollection> genParticles(new GenParticleCollection);

        // Run the generator and put the generated particles into the event.
        _primaryProtonGunGenerator->generate(*genParticles);
        evt.put(std::move(genParticles));

    }//produce()

}

DEFINE_ART_MODULE(mu2e::PrimaryProtonGun)
