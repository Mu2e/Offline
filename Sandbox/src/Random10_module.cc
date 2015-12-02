//
// Test multiple instances of engines in one module.
//
// Contact person Rob Kutschke
//

#include "SeedService/inc/SeedService.hh"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "CLHEP/Random/RandFlat.h"

#include <iostream>

using namespace std;

namespace mu2e {

  class Random10 : public art::EDAnalyzer {
  public:

    explicit Random10(fhicl::ParameterSet const& pset);

    virtual void analyze(art::Event const& e) override;

  private:

    art::RandomNumberGenerator::base_engine_t&     engine1_;
    art::RandomNumberGenerator::base_engine_t&     engine2_;
    CLHEP::RandFlat                                flat1_;
    CLHEP::RandFlat                                flat2_;

  };
}

mu2e::Random10::Random10(fhicl::ParameterSet const& pset):
  EDAnalyzer(pset),
  engine1_( createEngine(art::ServiceHandle<SeedService>()->getSeed("foo"),
                         "HepJamesRandom",
                         "foo") ),
  engine2_( createEngine(art::ServiceHandle<SeedService>()->getSeed("bar"),
                         "HepJamesRandom",
                         "bar") ),
  flat1_(engine1_),
  flat2_(engine2_){
  }

void mu2e::Random10::analyze( art::Event const& event ) {
  std::cout << "Event: "
            << event.id().event() << " "
            << flat1_.fire()  << " "
            << flat2_.fire()  << " "
            << std::endl;
}

DEFINE_ART_MODULE(mu2e::Random10);
