//
// Use TacerProduct to study emplace_back.
//
//
// Original author Rob Kutschke.
//

// Mu2e includes.
#include "Sandbox/inc/TracerProduct.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"

namespace mu2e {

  class MoveAwareTest00 : public art::EDAnalyzer {
  public:
    explicit MoveAwareTest00(fhicl::ParameterSet const& );

    void analyze( art::Event const& e);

  private:

  };

  MoveAwareTest00::MoveAwareTest00(fhicl::ParameterSet const& pset )
    : art::EDAnalyzer(pset){
  }

  void
  MoveAwareTest00::analyze(art::Event const& event) {

    const size_t n(3);
    std::vector<TracerProduct> v0;
    v0.reserve(n);

    mf::LogVerbatim("Tracing") << "Start first test: push_back as a baseline.";
    for ( size_t i=0; i<n; ++i ){
      v0.push_back( TracerProduct( i ) );
    }

    for ( size_t i=0; i< v0.size(); ++i ){
      mf::LogVerbatim("Tracing") << "Recap step 1: " << i << " " << v0[i];
    }
    v0.clear();
    v0.reserve(n);

    mf::LogVerbatim("Tracing") << "\n\nStart second test: check that emplace back works as advertised";
    for ( size_t i=0; i<n; ++i ){
      v0.emplace_back( i );
    }

    for ( size_t i=0; i< v0.size(); ++i ){
      mf::LogVerbatim("Tracing") << "Recap step 2: " << i << " " << v0[i];
    }

    mf::LogVerbatim("Tracing") << "\n\nStart third test: element by element copy from another vector, using push_back";
    std::vector<TracerProduct> v1;
    v1.reserve(n);
    for ( size_t i=0; i<n; ++i ){
      v1.push_back( v0[i] );
    }

    for ( size_t i=0; i< v1.size(); ++i ){
      mf::LogVerbatim("Tracing") << "Recap step 3: " << i << " " << v1[i];
    }

    mf::LogVerbatim("Tracing") << "\n\nStart fourth test: element by element copy from another vector, using emplace_back";
    std::vector<TracerProduct> v2;
    v2.reserve(n);
    for ( size_t i=0; i<n; ++i ){
      v2.emplace_back( v0[i] );
    }

    for ( size_t i=0; i< v2.size(); ++i ){
      mf::LogVerbatim("Tracing") << "Recap step 4: " << i << " " << v2[i];
    }

    mf::LogVerbatim("Tracing") << "\n\nStart fifth test: copy from another vector";
    std::vector<TracerProduct> v3(v0);
    for ( size_t i=0; i< v3.size(); ++i ){
      mf::LogVerbatim("Tracing") << "Recap step 5: " << i << " " << v3[i];
    }
    v0.clear();
    v1.clear();
    v2.clear();
    v3.clear();

    mf::LogVerbatim("Tracing") << "\n\nStart sixth test: see if never-used-again will trigger the move";
    TracerProduct a(42);
    mf::LogVerbatim("Tracing") << "a is: " << a;

    TracerProduct b = a;
    mf::LogVerbatim("Tracing") << "b is: " << b;

    TracerProduct c = std::move(a);
    mf::LogVerbatim("Tracing") << "c is: " << c;

  } // end of ::analyze.

}

DEFINE_ART_MODULE(mu2e::MoveAwareTest00)
