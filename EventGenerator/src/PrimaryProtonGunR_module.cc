
/*
  This is a Replicated Module.
  A plug-in for running PrimaryProtonGun-based event generator for running in MT art.
  It produces a GenParticleCollection of primary protons using the PrimaryProtonGun.
 
  These Collections are used in Mu2eG4_module.cc.
 
  Original author Lisa Goodenough

*/

// Mu2e includes.
#include "ConfigTools/inc/SimpleConfig.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"

// Particular generators that this code knows about.
#include "EventGenerator/inc/PrimaryProtonGun.hh"
#include "SeedService/inc/SeedService.hh"

// Includes from art and its toolchain.
#include "art/Framework/Core/ReplicatedProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes.
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace mu2e {

  class PrimaryProtonGunR : public art::ReplicatedProducer {

  public:

      explicit PrimaryProtonGunR(fhicl::ParameterSet const& pS, art::ProcessingFrame const& pF);
      // Accept compiler written d'tor.  Modules are never moved or copied.

      virtual void produce (art::Event& e, art::ProcessingFrame const& pF) override;
      virtual void beginRun(art::Run const& r, art::ProcessingFrame const& pF) override;

  private:

      // Name of the run-time configuration file.
      string _configfile;

      bool _allowReplacement;
      bool _messageOnReplacement;
      bool _messageOnDefault;
      int  _configStatsVerbosity;

      // Print final config file after all replacements.
      bool _printConfig;
    
      CLHEP::HepJamesRandom _engine;
      std::unique_ptr<PrimaryProtonGun> _primaryProtonGunGenerator;
      
      // Number of times BeginRun is called on this module
      int ncalls = 0;
      
  };

    PrimaryProtonGunR::PrimaryProtonGunR(fhicl::ParameterSet const& pSet, art::ProcessingFrame const& procFrame):
        art::ReplicatedProducer{pSet,procFrame},
        _configfile(            pSet.get<std::string>   ("inputfile")),
        _allowReplacement(      pSet.get<bool>          ("allowReplacement",      true)),
        _messageOnReplacement(  pSet.get<bool>          ("messageOnReplacement",  false)),
        _messageOnDefault(      pSet.get<bool>          ("messageOnDefault",      false)),
        _configStatsVerbosity(  pSet.get<int>           ("configStatsVerbosity",  0)),
        _printConfig(           pSet.get<bool>          ("printConfig",           false)),
        _engine{art::ServiceHandle<SeedService>{}->getSeed()}
    {
        produces<GenParticleCollection>();
    }
    
    void PrimaryProtonGunR::beginRun(art::Run const& run, art::ProcessingFrame const& procFrame){
        
    // The configuration of the PPG Generator does not change within a job.
    if ( ++ncalls > 1){
      mf::LogInfo("PrimaryProtonGunR")
        << "For Schedule: " << procFrame.scheduleID()
        << ", PrimaryProtonGunR Generator does not change state at beginRun.  Hope that's OK.";
      return;
    }

    // We don't want to print this out more than once,
    // regardless of the number of instances/schedules running.
    std::string schedID = std::to_string(procFrame.scheduleID().id());
        
    if ( schedID == "0"){
        cout << "Event generator configuration file: "
        << _configfile
        << "\n"
        << endl;
    }

    // Load the configuration, make modifications if required, and print if desired.
    SimpleConfig config(_configfile, _allowReplacement, _messageOnReplacement, _messageOnDefault );
    if ( _printConfig ){
      config.print(cout,"PrimaryProtonGunR: ");
    }
    config.printAllSummaries( cout, _configStatsVerbosity, "PrimaryProtonGunR: ");
     
        
    // Instantiate generator for this run.
    _primaryProtonGunGenerator = std::make_unique <PrimaryProtonGun>( _engine, run, config);
        
  }//beginRun

    void PrimaryProtonGunR::produce(art::Event& evt, art::ProcessingFrame const& procFrame) {
        
        // Make the collections to hold the output.
        unique_ptr<GenParticleCollection> genParticles(new GenParticleCollection);
        
        // Run the generator and put the generated particles into the event.
        _primaryProtonGunGenerator->generate(*genParticles);
        evt.put(std::move(genParticles));

    }//produce()

}

DEFINE_ART_MODULE(mu2e::PrimaryProtonGunR);
