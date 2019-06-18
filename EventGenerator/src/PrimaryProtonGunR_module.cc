
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

// Other external includes.
#include <boost/shared_ptr.hpp>

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

      virtual void produce (art::Event& e, art::ProcessingFrame const& pF);
      virtual void beginRun(art::Run const& r, art::ProcessingFrame const& pF);

  private:

      // Name of the run-time configuration file.
      string _configfile;

      // Control the behaviour of messages from the SimpleConfig object holding
      // the geometry parameters.
      bool _allowReplacement;
      bool _messageOnReplacement;
      bool _messageOnDefault;
      int  _configStatsVerbosity;

      // Print final config file after all replacements.
      bool _printConfig;
    
      //CLHEP::HepRandomEngine& _engine;
      CLHEP::HepJamesRandom _engine;
      
      // A collection of all of the generators that we will run.
      typedef boost::shared_ptr<GeneratorBase> GeneratorBasePtr;
      GeneratorBasePtr _primaryProtonGunGenerator;
      
      // Check for a valid configuration
      void checkConfig( const SimpleConfig& config );
      
      // Number of times BeginRun is called on this module
      int ncalls;
      
  };

    PrimaryProtonGunR::PrimaryProtonGunR(fhicl::ParameterSet const& pSet, art::ProcessingFrame const& procFrame):
        art::ReplicatedProducer{pSet,procFrame},
        _configfile(            pSet.get<std::string>   ("inputfile",             "generatorconfig.txt")),
        _allowReplacement(      pSet.get<bool>          ("allowReplacement",      true)),
        _messageOnReplacement(  pSet.get<bool>          ("messageOnReplacement",  false)),
        _messageOnDefault(      pSet.get<bool>          ("messageOnDefault",      false)),
        _configStatsVerbosity(  pSet.get<int>           ("configStatsVerbosity",  0)),
        _printConfig(           pSet.get<bool>          ("printConfig",           false)),
        _engine{art::ServiceHandle<SeedService>{}->getSeed()},
        ncalls(0)
    {
        produces<GenParticleCollection>();
    }


  // Run this at beginRun time in case any of them depend on geoemtry
  // information that may change with conditions information.
  // Otherwise this information could be computed in the c'tor.
    void PrimaryProtonGunR::beginRun(art::Run const& run, art::ProcessingFrame const& procFrame){
        
    //The configuration of the PPG Generator does not change within a job
    if ( ++ncalls > 1){
      mf::LogInfo("PrimaryProtonGunR")
        << "For Schedule: " << procFrame.scheduleID()
        << ", PrimaryProtonGunR Generator does not change state at beginRun.  Hope that's OK.";
      return;
    }

    //we don't want to print this out more than once, regardless of the number of instances
    static int instance(0);
    if ( instance == 0){
        cout << "Event generator configuration file: "
        << _configfile
        << "\n"
        << endl;
    }

    SimpleConfig config(_configfile, _allowReplacement, _messageOnReplacement, _messageOnDefault );
    checkConfig(config);

    if ( _printConfig ){
      config.print(cout,"PrimaryProtonGunR: ");
    }

    config.printAllSummaries( cout, _configStatsVerbosity, "PrimaryProtonGunR: ");
     
    //we need a different output data filename for each schedule
    //since BeginRun is called serially, we don't need to worry about a lock here
    //static int instance(0);
    // Instantiate generator for this run.
    _primaryProtonGunGenerator = GeneratorBasePtr( new PrimaryProtonGun( _engine, run, config, instance)) ;
    instance++;
        
  }//beginRun

    void PrimaryProtonGunR::produce(art::Event& evt, art::ProcessingFrame const& procFrame) {
        
        // Make the collections to hold the output.
        unique_ptr<GenParticleCollection> genParticles(new GenParticleCollection);
        
        // Run the generators.
        _primaryProtonGunGenerator->generate(*genParticles);
        // Put the generated particles into the event.
        evt.put(std::move(genParticles));

    }//produce()

    // Look for inconsistencies in the config file.
    void PrimaryProtonGunR::checkConfig( const SimpleConfig&  config){
      // There is nothing to do for this generator
    }//checkConfig


}


using mu2e::PrimaryProtonGunR;
DEFINE_ART_MODULE(PrimaryProtonGunR);
