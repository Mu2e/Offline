//
// Test the FooService
//
// $Id: FooTest_module.cc,v 1.1 2012/07/24 20:01:16 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/24 20:01:16 $
//
// Contact person Rob Kutschke
//

#include "Sandbox/inc/FooService.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

// C++ includes.
#include <iostream>

namespace mu2e {

  class FooTest : public art::EDAnalyzer {
  public:

    explicit FooTest(fhicl::ParameterSet const& pset);

    virtual void analyze(const art::Event& e);

  };

  FooTest::FooTest(fhicl::ParameterSet const& pset){
    std::cout << "FooTest: constructor:" << std::endl;
  }

  void FooTest::analyze( const art::Event& event ) {
    std::cout << "FooTest::analyze: " << event.id() << std::endl;
    art::ServiceHandle<FooService> foo;
    foo->poke();
  }
}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::FooTest);
