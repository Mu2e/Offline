//
//  The first example of a producer.
//
//  $Id: HelloProducer_plugin.cc,v 1.2 2011/03/04 23:34:25 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2011/03/04 23:34:25 $
//   
//  Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Framework includes.
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

// Mu2e includes.
#include "ToyDP/inc/ToyGenParticleCollection.hh"

using namespace std;

namespace mu2e {

  class HelloProducer : public edm::EDProducer {

  public:
    explicit HelloProducer(edm::ParameterSet const& pset){
      produces<ToyGenParticleCollection>();
    }

    void produce( edm::Event& event, edm::EventSetup const&);

  private:

  };

  void HelloProducer::produce( edm::Event& event, edm::EventSetup const&){

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
DEFINE_FWK_MODULE(HelloProducer);
