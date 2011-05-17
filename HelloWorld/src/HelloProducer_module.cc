//
//  The first example of a producer.
//
//  $Id: HelloProducer_module.cc,v 1.2 2011/05/17 22:22:46 wb Exp $
//  $Author: wb $
//  $Date: 2011/05/17 22:22:46 $
//   
//  Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"

// Mu2e includes.
#include "ToyDP/inc/ToyGenParticleCollection.hh"

using namespace std;

namespace mu2e {

  class HelloProducer : public art::EDProducer {

  public:
    explicit HelloProducer(fhicl::ParameterSet const& pset){
      produces<ToyGenParticleCollection>();
    }

    void produce( art::Event& event);

  private:

  };

  void HelloProducer::produce( art::Event& event){

    auto_ptr<ToyGenParticleCollection> genParticles(new ToyGenParticleCollection);

    CLHEP::Hep3Vector position(0.,0.,0.);
    CLHEP::HepLorentzVector momentum(50.,0.,0.,50.);
    double time(0.);

    genParticles->push_back( ToyGenParticle( PDGCode::gamma, 
                                             GenId::particleGun, 
                                             position,
                                             momentum,
                                             time)
                            );
    
    event.put( genParticles );

  }

} // end namespace mu2e

using mu2e::HelloProducer;
DEFINE_ART_MODULE(HelloProducer);
