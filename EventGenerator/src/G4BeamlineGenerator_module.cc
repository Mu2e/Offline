/*

  A plug_in for running G4Beamline-based event generator.

  Adds additional data collection to the event with extra data
  available in G4Beamline data file.

  The code is kept very similar to EventGenerator_plugin

*/

// Mu2e includes.
#include "ConfigTools/inc/SimpleConfig.hh"
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "SeedService/inc/SeedService.hh"

// Includes from art and its toolchain.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Principal/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Particular generators that this code knows about.
#include "EventGenerator/inc/FromG4BLFile.hh"

// Other external includes.
#include <boost/shared_ptr.hpp>

// C++ includes.
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace mu2e {

  class G4BeamlineGenerator : public art::EDProducer {

  public:

    explicit G4BeamlineGenerator(fhicl::ParameterSet const& pSet);
    // Accept compiler written d'tor.  Modules are never moved or copied.

    virtual void produce (art::Event& e);
    virtual void beginRun(art::Run&   r);

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

    // A collection of all of the generators that we will run.
    typedef  boost::shared_ptr<FromG4BLFile> GeneratorBasePtr;
    GeneratorBasePtr _generator;
    CLHEP::HepRandomEngine& _engine;

    // Check that configuration is internally consistent.
    void checkConfig( const SimpleConfig&  config);

  };

  G4BeamlineGenerator::G4BeamlineGenerator(fhicl::ParameterSet const& pSet):
    EDProducer{pSet},
    _configfile(           pSet.get<std::string>("inputfile",            "generatorconfig.txt")),
    _allowReplacement(     pSet.get<bool>       ("allowReplacement",     true)),
    _messageOnReplacement( pSet.get<bool>       ("messageOnReplacement", true)),
    _messageOnDefault(     pSet.get<bool>       ("messageOnDefault",      false)),
    _configStatsVerbosity( pSet.get<int>        ("configStatsVerbosity",  0)),
    _printConfig(          pSet.get<bool>       ("printConfig",           false)),
    // Provide a common engine for the generators to use via the service
    _engine{createEngine(art::ServiceHandle<SeedService>{}->getSeed())}
  {
    produces<GenParticleCollection>();
    produces<G4BeamlineInfoCollection>();
  }


  // At beginRun time, update any derived geometry information.
  void G4BeamlineGenerator::beginRun( art::Run &run){

    static int ncalls(0);
    if ( ++ncalls > 1){
      mf::LogInfo("G4BeamlineGenerator")
        << "G4BeamlineGenerator does not change state at beginRun.  Hope that's OK.";
      return;
    }

    cout << "Event generator configuration file: "
         << _configfile
         << "\n"
         << endl;

    SimpleConfig config(_configfile, _allowReplacement, _messageOnReplacement, _messageOnDefault );
    checkConfig(config);

    if ( _printConfig ){
      config.print(cout, "G4blGen: " );
    }

    // Instantiate generators for this run.
    _generator = GeneratorBasePtr(new FromG4BLFile{_engine, run, config});

    config.printAllSummaries( cout, _configStatsVerbosity, "G4blGen: ");
  }

  void
  G4BeamlineGenerator::produce(art::Event& evt) {

    // Make the collection to hold the output.
    unique_ptr<GenParticleCollection> genParticles(new GenParticleCollection);
    unique_ptr<G4BeamlineInfoCollection> extraData(new G4BeamlineInfoCollection);

    _generator->generate(*genParticles,&(*extraData));

    // Put the generated particles into the event.
    evt.put(std::move(genParticles));
    evt.put(std::move(extraData));

  }

  // Look for inconsistencies in the config file.
  void G4BeamlineGenerator::checkConfig( const SimpleConfig&  config){

    // There is nothing to do for this generator

  }


}


using mu2e::G4BeamlineGenerator;
DEFINE_ART_MODULE(G4BeamlineGenerator);
