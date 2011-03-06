
/*

  A plug_in for running G4Beamline-based event generator.

  Adds additional data collection to the event with extra data 
  available in G4Beamline data file. 

  The code is kept very similar to EventGenerator_plugin

*/

// C++ includes.
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

// Framework includes.
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// Mu2e includes.
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "ToyDP/inc/GenId.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/G4BeamlineInfoCollection.hh"

// Particular generators that this code knows about.
#include "EventGenerator/inc/FromG4BLFile.hh"

// Other external includes.
#include <boost/shared_ptr.hpp>

using namespace std;

namespace mu2e {

  class G4BeamlineGenerator : public edm::EDProducer {

  public:

    explicit G4BeamlineGenerator(edm::ParameterSet const& pSet):
      _configfile(pSet.getUntrackedParameter<std::string>("inputfile","generatorconfig.txt"))
    {
      // A placeholder until I make a real data product.
      produces<ToyGenParticleCollection>();
      produces<G4BeamlineInfoCollection>();

      // Print generators for which Id's are defined.
      //GenId::printAll();

      // Provide a common engine for the generators to use via the service
      createEngine( get_seed_value(pSet) );
    }

    virtual ~G4BeamlineGenerator() { }

    virtual void produce(edm::Event& e, edm::EventSetup const& c);

    virtual void beginRun(edm::Run &r, edm::EventSetup const& eSetup );

    static void fillDescription(edm::ParameterSetDescription& iDesc,
                                string const& moduleLabel) {
      iDesc.setAllowAnything();
    }

  private:

    // Name of the run-time configuration file.
    string _configfile;

    // A collection of all of the generators that we will run.
    typedef  boost::shared_ptr<FromG4BLFile> GeneratorBasePtr;
    GeneratorBasePtr _generator;

    void checkConfig( const SimpleConfig&  config);

  };

  // At beginRun time, update any derived geometry information.
  void G4BeamlineGenerator::beginRun( edm::Run &run, edm::EventSetup const& eSetup ){

    static int ncalls(0);
    if ( ++ncalls > 1){
      edm::LogInfo("G4BeamlineGenerator")
        << "G4BeamlineGenerator does not change state at beginRun.  Hope that's OK.";
      return;
    }

    edm::LogInfo log("G4BeamlineGenerator");
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
  G4BeamlineGenerator::produce(edm::Event& evt, edm::EventSetup const&) {

    // Make the collection to hold the output.
    auto_ptr<ToyGenParticleCollection> genParticles(new ToyGenParticleCollection);
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
DEFINE_FWK_MODULE(G4BeamlineGenerator);
