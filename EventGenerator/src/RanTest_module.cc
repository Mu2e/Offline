/*

 Plugin to test access of Random number generator

*/

// Mu2e includes.
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "SeedService/inc/SeedService.hh"

// Includes from art and its tool chain
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "fhiclcpp/ParameterSet.h"

// C++ includes.
#include <iostream>

using namespace std;

namespace mu2e {

  class RanTest : public art::EDProducer {

  public:

    explicit RanTest(fhicl::ParameterSet const& pSet);
    virtual ~RanTest();

    virtual void produce(art::Event& e);

  private:

  };

  RanTest::RanTest(fhicl::ParameterSet const& pSet) : EDProducer{pSet} {

    produces<GenParticleCollection>();

    // Provide a common engine for the generators to use via the service
    createEngine( art::ServiceHandle<SeedService>()->getSeed() );
  }
  RanTest::~RanTest() { }


  void
  RanTest::produce(art::Event& evt) {

    cout << "Creating an empty collection ... " << endl;
    unique_ptr<GenParticleCollection> genParticles(new GenParticleCollection);
    evt.put(std::move(genParticles));

  }
}


using mu2e::RanTest;
DEFINE_ART_MODULE(RanTest);
