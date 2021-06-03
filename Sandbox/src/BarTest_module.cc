//
// Test the BarService
//
//
// Contact person Rob Kutschke
//

#include "Sandbox/inc/BarService.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// C++ includes.
#include <iostream>

namespace mu2e {

  class BarTest : public art::EDAnalyzer {
  public:

    explicit BarTest(fhicl::ParameterSet const& pset);

    virtual void analyze(const art::Event& e);

  };

  BarTest::BarTest(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset){
    std::cout << "BarTest: constructor:" << std::endl;
  }

  void BarTest::analyze( const art::Event& event ) {
    std::cout << "BarTest::analyze: " << event.id() << std::endl;
    art::ServiceHandle<BarService> foo;
    foo->poke();
  }
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::BarTest);
