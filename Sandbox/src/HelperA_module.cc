//
// Explore the One Definition Rule
//
// Contact person Rob Kutschke
//

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

// C++ includes.
#include <iostream>

namespace mu2e {

  class Helper{
  public:
    void doit( double a );
  private:
  };

  class HelperA : public art::EDAnalyzer {
  public:

    explicit HelperA(fhicl::ParameterSet const& pset);

    virtual void analyze(const art::Event& e);

    Helper h_;

  };

  HelperA::HelperA(fhicl::ParameterSet const& pset):
    art::EDAnalyzer(pset),
    h_(){
  }

  void HelperA::analyze( const art::Event& ) {
    h_.doit(2.);
  }

  void Helper::doit( double a ){
    std::cout << "This is the class from HelperA: " << a << std::endl;
  }

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::HelperA);
