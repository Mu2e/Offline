
/*
  This is a Replicated Module.
  A plug-in for running PrimaryProtonGun-based event generator for running in MT art.
  It produces a GenParticleCollection of primary protons using the PrimaryProtonGun.
 
  These Collections are used in Mu2eG4_module.cc.
 
 Original author: Lisa Goodenough
  Date: 2020/02/19

*/

// Mu2e includes.
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

// Particular generators that this code knows about.
#include "EventGenerator/inc/PrimaryProtonGunImpl.hh"
#include "SeedService/inc/SeedService.hh"

// Includes from art and its toolchain.
#include "art/Framework/Core/ReplicatedProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes.
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

namespace mu2e {

  class PrimaryProtonGun : public art::ReplicatedProducer {

  public:
      
      typedef art::ReplicatedProducer::Table<PrimaryProtonGunImpl::PrimaryProtonGunConfig> Parameters;
      
      explicit PrimaryProtonGun(Parameters const& params, art::ProcessingFrame const& pF);
      // Accept compiler written d'tor.  Modules are never moved or copied.

      virtual void produce(art::Event& e, art::ProcessingFrame const& pF) override;
      virtual void beginRun(art::Run const& r, art::ProcessingFrame const& pF) override;

  private:

      const PrimaryProtonGunImpl::PrimaryProtonGunConfig _config;

      //This engine implementation is only necessary for art versions below v3_02_06
      //the fix in v3_02_06 allows The RandomNumberService to be used in a Replicated Module
      //CLHEP::HepJamesRandom _engine;
      CLHEP::HepRandomEngine& _engine;
      std::unique_ptr<PrimaryProtonGunImpl> _primaryProtonGunGenerator;
      
      // Number of times BeginRun is called on this module
      int ncalls = 0;
      
  };

    PrimaryProtonGun::PrimaryProtonGun(Parameters const& params, art::ProcessingFrame const& procFrame):
        art::ReplicatedProducer{params,procFrame},
        _config(params()),
        _engine{createEngine(art::ServiceHandle<SeedService>{}->getSeed())}
    {
        produces<GenParticleCollection>();
    }
    
    void PrimaryProtonGun::beginRun(art::Run const& run, art::ProcessingFrame const& procFrame){
        
    // The configuration of the PPG Generator does not change within a job.
    if ( ++ncalls > 1){
      mf::LogInfo("PrimaryProtonGun")
        << "For Schedule: " << procFrame.scheduleID()
        << ", PrimaryProtonGun Generator does not change state at beginRun.  Hope that's OK.";
      return;
    }

    // Instantiate generator for this run.
    _primaryProtonGunGenerator = std::make_unique <PrimaryProtonGunImpl>( _engine, _config);
        
  }//beginRun

    void PrimaryProtonGun::produce(art::Event& evt, art::ProcessingFrame const& procFrame) {
        
        // Make the collections to hold the output.
        unique_ptr<GenParticleCollection> genParticles(new GenParticleCollection);
        
        // Run the generator and put the generated particles into the event.
        _primaryProtonGunGenerator->generate(*genParticles);
        evt.put(std::move(genParticles));

    }//produce()

}

DEFINE_ART_MODULE(mu2e::PrimaryProtonGun);
