
/*

  A plug-in for running PrimaryProtonGun-based event generator for running in MT mode.
  It produces a collection of GenParticleCollections of
  primary protons using the PrimaryProtonGun.  These Collections are
  used in the stashes in MT mode of Mu2eG4_module.cc.
 
  The code is tightly based on the EventGenerator_plugin

  Original author Lisa Goodenough

*/

// Mu2e includes.
#include "ConfigTools/inc/SimpleConfig.hh"
#include "MCDataProducts/inc/GenId.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/GenParticleCollections.hh"

// Particular generators that this code knows about.
#include "EventGenerator/inc/PrimaryProtonGun.hh"
#include "SeedService/inc/SeedService.hh"

// Includes from art and its toolchain.
#include "art/Framework/Core/EDProducer.h"
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

  class PrimaryProtonGunN : public art::EDProducer {

  public:

      explicit PrimaryProtonGunN(fhicl::ParameterSet const& pSet);
      // Accept compiler written d'tor.  Modules are never moved or copied.

      virtual void produce (art::Event& e);
      virtual void beginRun(art::Run&   r);
      void resetStash();
      void fillStash();

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
      
      CLHEP::HepRandomEngine& _engine;

      // A collection of all of the generators that we will run.
      typedef boost::shared_ptr<GeneratorBase> GeneratorBasePtr;
      GeneratorBasePtr _primaryProtonGunGenerator;
      
      size_t stashSize_;
      size_t finger_;
      GenParticleCollections stash_;
      
      // Check for a valid configuration
      void checkConfig( const SimpleConfig& config );

  };

  PrimaryProtonGunN::PrimaryProtonGunN(fhicl::ParameterSet const& pSet):
        art::EDProducer{pSet},
        _configfile(            pSet.get<std::string>   ("inputfile",             "generatorconfig.txt")),
        _allowReplacement(      pSet.get<bool>          ("allowReplacement",      true)),
        _messageOnReplacement(  pSet.get<bool>          ("messageOnReplacement",  false)),
        _messageOnDefault(      pSet.get<bool>          ("messageOnDefault",      false)),
        _configStatsVerbosity(  pSet.get<int>           ("configStatsVerbosity",  0)),
        _printConfig(           pSet.get<bool>          ("printConfig",           false)),
        // A common random engine for the generator to use.
        _engine{createEngine(art::ServiceHandle<SeedService>{}->getSeed())},
        stashSize_(             pSet.get<size_t>        ("stashSize",             1)),
        finger_(stashSize_),
        stash_(stashSize_)    
    {
        produces<GenParticleCollection>();
        produces<GenParticleCollections>();
    }


  // Run this at beginRun time in case any of them depend on geoemtry
  // information that may change with conditions information.
  // Otherwise this information could be computed in the c'tor.
  void PrimaryProtonGunN::beginRun( art::Run &run){

    static int ncalls(0);
    if ( ++ncalls > 1){
      mf::LogInfo("PrimaryProtonGunN")
        << "PrimaryProtonGunN Generator does not change state at beginRun.  Hope that's OK.";
      return;
    }

    cout << "Event generator configuration file: "
         << _configfile
         << "\n"
         << endl;

    SimpleConfig config(_configfile, _allowReplacement, _messageOnReplacement, _messageOnDefault );
    checkConfig(config);

    if ( _printConfig ){
      config.print(cout,"PrimaryProtonGunGenN: ");
    }

    // Instantiate generator for this run.
    _primaryProtonGunGenerator = GeneratorBasePtr( new PrimaryProtonGun( _engine, run, config )) ;
 
    config.printAllSummaries( cout, _configStatsVerbosity, "PrimaryProtonGunGenN: ");

  }//beginRun

    void PrimaryProtonGunN::produce(art::Event& evt) {

        // Make the collections to hold the output.
        //unique_ptr<GenParticleCollection> genParticles(new GenParticleCollection);
        // Run the generators.
        //_primaryProtonGunGenerator->generate(*genParticles);
        // Put the generated particles into the event.
        //evt.put(std::move(genParticles));
        
        //*****************************************************************
        // On event 0, N, 2N, ... fill the stash with N events
        // and place a copy of the stash in the event.
        if ( finger_ == stash_.size() ){
            resetStash();
            fillStash();
            finger_=0;
            auto allEvents = std::make_unique<GenParticleCollections>(stash_);
            evt.put(std::move(allEvents));
        } else{
            auto allEvents = std::make_unique<GenParticleCollections>();
            evt.put(std::move(allEvents));
        }
        
        // On every event, copy the next GenParticleCollection out of the stash and
        // put it into the event.
        auto gens = std::make_unique<GenParticleCollection>(stash_[finger_++]);
        evt.put(std::move(gens));

    }//produce()

    
    void PrimaryProtonGunN::resetStash() {
        for ( size_t i=0; i<stashSize_; ++i){
            stash_[i].clear();
        }
    }
    
    void PrimaryProtonGunN::fillStash() {
        
        for ( size_t i=0; i<stashSize_; ++i){
            _primaryProtonGunGenerator->generate(stash_[i]);
        }
    }

    // Look for inconsistencies in the config file.
    void PrimaryProtonGunN::checkConfig( const SimpleConfig&  config){
      // There is nothing to do for this generator
    }


}


using mu2e::PrimaryProtonGunN;
DEFINE_ART_MODULE(PrimaryProtonGunN);
