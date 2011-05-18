
/*

  A plug_in for running G4Beamline-based event generator.

  Adds additional data collection to the event with extra data
  available in G4Beamline data file.

  The code is kept very similar to EventGenerator_plugin

*/

// C++ includes.
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/Handle.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "ToyDP/inc/G4BeamlineInfoCollection.hh"
#include "ToyDP/inc/GenId.hh"
#include "ToyDP/inc/GenParticleCollection.hh"

// Particular generators that this code knows about.
#include "EventGenerator/inc/FromG4BLFile.hh"

// Other external includes.
#include <boost/shared_ptr.hpp>

using namespace std;

namespace mu2e {

  class G4BeamlineGenerator : public art::EDProducer {

  public:

    explicit G4BeamlineGenerator(fhicl::ParameterSet const& pSet):
      _configfile(pSet.get<std::string>("inputfile","generatorconfig.txt"))
    {
      // A placeholder until I make a real data product.
      produces<GenParticleCollection>();
      produces<G4BeamlineInfoCollection>();

      // Print generators for which Id's are defined.
      //GenId::printAll();

      // Provide a common engine for the generators to use via the service
      createEngine( get_seed_value(pSet) );
    }

    virtual ~G4BeamlineGenerator() { }

    virtual void produce(art::Event& e);

    virtual void beginRun(art::Run &r);

  private:

    // Name of the run-time configuration file.
    string _configfile;

    // A collection of all of the generators that we will run.
    typedef  boost::shared_ptr<FromG4BLFile> GeneratorBasePtr;
    GeneratorBasePtr _generator;

    void checkConfig( const SimpleConfig&  config);

  };

  // At beginRun time, update any derived geometry information.
  void G4BeamlineGenerator::beginRun( art::Run &run){

    static int ncalls(0);
    if ( ++ncalls > 1){
      mf::LogInfo("G4BeamlineGenerator")
        << "G4BeamlineGenerator does not change state at beginRun.  Hope that's OK.";
      return;
    }

    mf::LogInfo log("G4BeamlineGenerator");
    log << "Event generator configuration file: "
        << _configfile
        << "\n\n";

    SimpleConfig config(_configfile);
    checkConfig(config);

    if ( config.getBool("printConfig",false) ){
      log << config;
    }

    if ( config.getBool("printConfigStats",false) ){
      // Work around absence of << operator for this print method.
      ostringstream os;
      config.printStatistics(os);
      log << os.str();
    }

    // Instantiate generators for this run.
    _generator = GeneratorBasePtr(new FromG4BLFile(run, config));

  }

  void
  G4BeamlineGenerator::produce(art::Event& evt) {

    // Make the collection to hold the output.
    auto_ptr<GenParticleCollection> genParticles(new GenParticleCollection);
    auto_ptr<G4BeamlineInfoCollection> extraData(new G4BeamlineInfoCollection);

    _generator->generate(*genParticles,&(*extraData));

    // Put the generated particles into the event.
    evt.put(genParticles);
    evt.put(extraData);

  }

  // Look for inconsistencies in the config file.
  void G4BeamlineGenerator::checkConfig( const SimpleConfig&  config){

    // There is nothing to do for this generator

  }


}


using mu2e::G4BeamlineGenerator;
DEFINE_ART_MODULE(G4BeamlineGenerator);
