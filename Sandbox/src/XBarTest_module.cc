//
// Test the XBarService
//
//
// Contact person Rob Kutschke
//

#include "Sandbox/inc/XBarService.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// C++ includes.
#include <iostream>

namespace mu2e {

  class XBarTest : public art::EDAnalyzer {
  public:

    explicit XBarTest(fhicl::ParameterSet const& pset);

    virtual void analyze(const art::Event& e);

  };

  XBarTest::XBarTest(fhicl::ParameterSet const& pset)
    : art::EDAnalyzer(pset){
    std::cout << "XBarTest: constructor:" << std::endl;
  }

  void XBarTest::analyze( const art::Event& event ) {
    std::cout << "XBarTest::analyze: " << event.id() << std::endl;
    art::ServiceHandle<XBarService> foo;
    foo->poke();
  }
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::XBarTest);
