/*

 Plugin to test access of Random number generator

*/

// Mu2e includes.
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/SpectrumConfig.hh"
#include "Offline/SeedService/inc/SeedService.hh"

// Includes from art and its tool chain
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
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
    virtual void endSubRun(art::SubRun& sr) override;

  private:

  };

  RanTest::RanTest(fhicl::ParameterSet const& pSet) : EDProducer{pSet} {

    produces<GenParticleCollection>();
    produces<SpectrumConfig, art::InSubRun>();

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

  void RanTest::endSubRun(art::SubRun& sr) {
    auto config = std::make_unique<SpectrumConfig>();
    config->type_ = SpectrumConfig::Type::kOther;
    sr.put(std::move(config), art::fullSubRun());
  }
}


using mu2e::RanTest;
DEFINE_ART_MODULE(RanTest)
