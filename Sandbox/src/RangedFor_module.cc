//
// Use the TrackerProduct to check the behaviour of a ranged for.
//
//
// Original author Rob Kutschke
//

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

// Mu2e includes.
#include "Sandbox/inc/TracerProduct.hh"

#include <iostream>

using namespace std;

namespace mu2e {

  class RangedFor : public art::EDAnalyzer {
  public:

    explicit RangedFor(fhicl::ParameterSet const& pset);

    void analyze( art::Event const& e) override;

  private:

  };

  RangedFor::RangedFor(fhicl::ParameterSet const& pset ):
    art::EDAnalyzer(pset){
  }

  void RangedFor::analyze(art::Event const& event) {

    std::vector<TracerProduct> v;
    v.reserve(3);

    v.emplace_back( 11 );
    v.emplace_back( 22 );
    v.emplace_back( 33 );

    // This does not make extra copies.
    for ( auto const& t : v ){
      cerr << "Ref: " << t << endl;
    }

    // This does make extra copies.
    for ( auto  t : v ){
      cerr << "NonRef: " << t << endl;
    }


  } // end RangedFor::analyze

}  // end namespace mu2e

DEFINE_ART_MODULE(mu2e::RangedFor)
