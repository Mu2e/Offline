//
//  The first example of a producer.
//
//
//  Original author Rob Kutschke
//

// C++ includes.
#include <iostream>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

// Mu2e includes.
#include "MCDataProducts/inc/GenParticleCollection.hh"

using namespace std;

namespace mu2e {

  class HelloProducer : public art::EDProducer {

  public:
    explicit HelloProducer(fhicl::ParameterSet const& pset) :
      art::EDProducer{pset}
    {
      produces<GenParticleCollection>();
    }

    void produce( art::Event& event);

  private:

  };

  void HelloProducer::produce( art::Event& event){

    std::unique_ptr<GenParticleCollection> genParticles(new GenParticleCollection);

    CLHEP::Hep3Vector position(0.,0.,0.);
    CLHEP::HepLorentzVector momentum(50.,0.,0.,50.);
    double time(0.);

    genParticles->push_back( GenParticle( PDGCode::gamma,
                                             GenId::particleGun,
                                             position,
                                             momentum,
                                             time)
                            );

    event.put( std::move(genParticles) );

  }

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::HelloProducer);
