
/*

  A plug-in for running PrimaryProtonGun-based event generator for running in MT mode.
  It produces a GenParticleCollection of primary protons using the PrimaryProtonGun.
 
  These Collections are used in Mu2eG4_module.cc.
 
  The code is tightly based on the EventGenerator_plugin

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

  };

    PrimaryProtonGunR::PrimaryProtonGunR(fhicl::ParameterSet const& pSet, art::ProcessingFrame const& procFrame):
        art::ReplicatedProducer{pSet,procFrame},
        _configfile(            pSet.get<std::string>   ("inputfile",             "generatorconfig.txt")),
        _allowReplacement(      pSet.get<bool>          ("allowReplacement",      true)),
        _messageOnReplacement(  pSet.get<bool>          ("messageOnReplacement",  false)),
        _messageOnDefault(      pSet.get<bool>          ("messageOnDefault",      false)),
        _configStatsVerbosity(  pSet.get<int>           ("configStatsVerbosity",  0)),
        _printConfig(           pSet.get<bool>          ("printConfig",           false)),
        // A common random engine for the generator to use.
        //_engine{createEngine(procFrame.serviceHandle<SeedService>{}->getSeed())}
        _engine{art::ServiceHandle<SeedService>{}->getSeed()}
        //_engine{createEngine(art::ServiceHandle<SeedService>{}->getSeed())}
    {
        produces<GenParticleCollection>();
        std::cout << "Constructing PrimaryProtonGunR #" << procFrame.scheduleID() << std::endl;
    }


  // Run this at beginRun time in case any of them depend on geoemtry
  // information that may change with conditions information.
  // Otherwise this information could be computed in the c'tor.
    void PrimaryProtonGunR::beginRun(art::Run const& run, art::ProcessingFrame const& procFrame){

    static int ncalls(0);
    if ( ++ncalls > 1){
      mf::LogInfo("PrimaryProtonGunR")
        << "PrimaryProtonGunR Generator does not change state at beginRun.  Hope that's OK.";
      return;
    }

    cout << "Event generator configuration file: "
         << _configfile
         << "\n"
         << endl;

    SimpleConfig config(_configfile, _allowReplacement, _messageOnReplacement, _messageOnDefault );
    checkConfig(config);

    if ( _printConfig ){
      config.print(cout,"PrimaryProtonGunR: ");
    }

    // Instantiate generator for this run.
    _primaryProtonGunGenerator = GeneratorBasePtr( new PrimaryProtonGun( _engine, run, config )) ;
 
    config.printAllSummaries( cout, _configStatsVerbosity, "PrimaryProtonGunR: ");

  }//beginRun

    void PrimaryProtonGunR::produce(art::Event& evt, art::ProcessingFrame const& procFrame) {

        std::cout << "THIS IS Schedule " << procFrame.scheduleID() << " at event " << evt.id() << std::endl;

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
    }


}


using mu2e::PrimaryProtonGunR;
DEFINE_ART_MODULE(PrimaryProtonGunR);
