
/*

  A plug_in for running a variety of event generators.

  $Id: EventGenerator_plugin.cc,v 1.3 2009/12/09 18:55:18 rhbob Exp $
  $Author: rhbob $
  $Date: 2009/12/09 18:55:18 $

  Original author Rob Kutschke

  Eventually this will support a variety of event generators, controllable
  from the run time configuration.  A given call might invoke one or more
  of these generators.
  
  1) A full featured single particle gun.
  2) Single conversion track, uniformly from the targets.
  3) (Emax-E)**5 DIO model.
  4) Other DIO models.
  5) protons, neutrons, gammas and nuclear fragments from muon capture.
  6) Mockups of pion capture on nuclei and of pion and muon decay in flight.
  I say mock-ups because I see this starting from an known pion and muon
  flux distributions, not by starting from a pion or a muon entering
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
#include "ToyDP/inc/ToyGenParticleCollection.hh"
#include "ToyDP/inc/GenId.hh"

// Particular generators that this code knows about.
#include "EventGenerator/inc/ParticleGun.hh"
#include "EventGenerator/inc/ConversionGun.hh"
#include "EventGenerator/inc/CosmicToy.hh"
#include "EventGenerator/inc/PiCapture.hh"
#include "EventGenerator/inc/DecayInOrbitGun.hh"

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

  };

  // At beginRun time, update any derived geometry information.
  void EventGenerator::beginRun( edm::Run &run, edm::EventSetup const& eSetup ){
    
    SimpleConfig config(_configfile);

    // Change this to modify rather than delete and make an new one??
    // Also only instantiate ones that are present in the config file.

    // Delete generators from the previous run.
    _generators.clear();

    // Instantiate generators for this run.
    _generators.push_back( GeneratorBasePtr( new ConversionGun(       run, config)) );
    _generators.push_back( GeneratorBasePtr( new CosmicToy(           run, config)) );
    _generators.push_back( GeneratorBasePtr( new PiCapture(           run, config)) );
    _generators.push_back( GeneratorBasePtr( new DecayInOrbitGun(     run, config)) );
    
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


}


using mu2e::EventGenerator;
DEFINE_FWK_MODULE(EventGenerator);
