
/*

  A plug_in for running a variety of event generators.

  $Id: EventGenerator_plugin.cc,v 1.11 2010/05/17 21:47:33 genser Exp $
  $Author: genser $
  $Date: 2010/05/17 21:47:33 $

  Original author Rob Kutschke

  Eventually this will support a variety of event generators, controllable
  from the run time configuration.  A given call might invoke one or more
  of these generators.
  
  1) A full featured single particle gun.
  2) Single conversion track, uniformly from the targets.
  3) (Emax-E)**5 DIO model.
  4) Other DIO models.
  5) protons, neutrons, gammas and nuclear fragments from muon capture.
  6) Mockups of CLHEP::pion capture on nuclei and of CLHEP::pion and muon decay in flight.
  I say mock-ups because I see this starting from an known CLHEP::pion and muon
  flux distributions, not by starting from a CLHEP::pion or a muon entering
  the DS.
  7) Simplified models of cosmics.
  
  At present I expect that the highest fidelity generation of cosmics will be
  done by running an external generator and then reading "events" from the output
  of that generator.  Perhaps the merge will be done in this module, perhaps
  it will be done in a separate module?

*/

// C++ includes.
#include <iostream>
#include <cassert>
#include <stdexcept>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>

// Framework includes.
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e includes.
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eUtilities/inc/requireUniqueKey.hh"
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/GenId.hh"

// Particular generators that this code knows about.
#include "EventGenerator/inc/ParticleGun.hh"
#include "EventGenerator/inc/ConversionGun.hh"
#include "EventGenerator/inc/CosmicToy.hh"
#include "EventGenerator/inc/CosmicDYB.hh"
#include "EventGenerator/inc/PiCapture.hh"
#include "EventGenerator/inc/DecayInOrbitGun.hh"
#include "EventGenerator/inc/EjectedProtonGun.hh"
#include "EventGenerator/inc/PiEplusNuGun.hh"
#include "EventGenerator/inc/PrimaryProtonGun.hh"

// Other external includes.
#include <boost/shared_ptr.hpp>

using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class EventGenerator : public edm::EDProducer {
    
  public:

    explicit EventGenerator(edm::ParameterSet const& pSet):
      _configfile(pSet.getUntrackedParameter<std::string>("inputfile","generatorconfig.txt"))
    {
      // A placeholder until I make a real data product.
      produces<ToyGenParticleCollection>();

      // Print generators for which Id's are defined.
      //GenId::printAll();

    }

    virtual ~EventGenerator() { }

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
    typedef  boost::shared_ptr<GeneratorBase> GeneratorBasePtr;    
    std::vector<GeneratorBasePtr> _generators;

    void checkConfig( const SimpleConfig&  config);

  };

  // At beginRun time, update any derived geometry information.
  void EventGenerator::beginRun( edm::Run &run, edm::EventSetup const& eSetup ){

    static int ncalls(0);
    if ( ++ncalls > 1){
      edm::LogInfo("EventGenerator") 
	<< "EventGenerator does not change state at beginRun.  Hope that's OK.";
    }
    
    SimpleConfig config(_configfile);
    checkConfig(config);
    
    // Change this to modify rather than delete and make an new one??

    // Delete generators from the previous run.
    _generators.clear();

    // Which generators will we run?
    bool doConv                 = config.getBool( "conversionGun.do", 1);
    bool doParticleGun          = config.getBool( "particleGun.do",   0);
    bool doCosmicToy            = config.getBool( "cosmictoy.do", 0);
    bool doCosmicDYB            = config.getBool( "cosmicDYB.do", 0);
    bool doPiCapture            = config.getBool( "picapture.do", 0);
    bool doEjectedProton        = config.getBool( "ejectedProtonGun.do", 0);
    bool doDIO                  = config.getBool( "decayinorbitGun.do", 0);
    bool doPiEplusNu            = config.getBool( "piEplusNuGun.do", 0);
    bool doPrimaryProtonGun     = config.getBool( "primaryProtonGun.do", 0);

    // Instantiate generators for this run.
    if ( doParticleGun)          _generators.push_back( GeneratorBasePtr( new ParticleGun(      run, config)) );
    if ( doConv)                 _generators.push_back( GeneratorBasePtr( new ConversionGun(    run, config)) );
    if ( doCosmicToy)            _generators.push_back( GeneratorBasePtr( new CosmicToy(        run, config)) );
    if ( doCosmicDYB)            _generators.push_back( GeneratorBasePtr( new CosmicDYB(        run, config)) );
    if ( doPiCapture)            _generators.push_back( GeneratorBasePtr( new PiCapture(        run, config)) );
    if ( doDIO)                  _generators.push_back( GeneratorBasePtr( new DecayInOrbitGun(  run, config)) );
    if ( doEjectedProton)        _generators.push_back( GeneratorBasePtr( new EjectedProtonGun( run, config)) );
    if ( doPiEplusNu)            _generators.push_back( GeneratorBasePtr( new PiEplusNuGun(     run, config)) );
    if ( doPrimaryProtonGun)     _generators.push_back( GeneratorBasePtr( new PrimaryProtonGun(     run, config)) );

    if ( _generators.size() == 0 ){
      edm::LogWarning("CONTROL")
        << "EventGenerator has no generators enabled. Hope that's OK.";
    }

  }
  
  void
  EventGenerator::produce(edm::Event& evt, edm::EventSetup const&) {

    // Make the collection to hold the output.
    auto_ptr<ToyGenParticleCollection> genParticles(new ToyGenParticleCollection);

    // Run all of the registered generators.
    for ( std::vector<GeneratorBasePtr>::const_iterator 
            i=_generators.begin(),
            e=_generators.end(); 
          i !=e; ++i ){
      (*i)->generate(*genParticles);
    }

    // Put the generated particles into the event.
    edm::OrphanHandle<ToyGenParticleCollection> q = evt.put(genParticles);

  }

  // Look for inconsistencies in the config file.
  void EventGenerator::checkConfig( const SimpleConfig&  config){

    // The known cosmic ray generators.
    vector<string> keys;
    if ( keys.size() == 0 ){
      keys.push_back("cosmictoy.do");
      keys.push_back("cosmicDYB.do");
      keys.push_back("cosmic.do");
    }

    // Require that 0 or 1 of the generators to be present.
    requireUniqueKey( keys, config );

  }


}


using mu2e::EventGenerator;
DEFINE_FWK_MODULE(EventGenerator);
