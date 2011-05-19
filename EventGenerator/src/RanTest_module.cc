
/*

 Plugin to test access of Random number generator

*/

// C++ includes.
#include <iostream>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "fhiclcpp/ParameterSet.h"

// Mu2e includes.
#include "ToyDP/inc/GenParticleCollection.hh"

using namespace std;

namespace mu2e {

  class RanTest : public art::EDProducer {

  public:

    explicit RanTest(fhicl::ParameterSet const& pSet);
    virtual ~RanTest();

    virtual void produce(art::Event& e);

  private:

  };

  RanTest::RanTest(fhicl::ParameterSet const& pSet){

    produces<GenParticleCollection>();

    // Provide a common engine for the generators to use via the service
    createEngine( get_seed_value(pSet) );
  }
  RanTest::~RanTest() { }


  void
  RanTest::produce(art::Event& evt) {

    cout << "Creating an empty collection ... " << endl;
    auto_ptr<GenParticleCollection> genParticles(new GenParticleCollection);
    evt.put(genParticles);

  }
}


using mu2e::RanTest;
DEFINE_ART_MODULE(RanTest);
