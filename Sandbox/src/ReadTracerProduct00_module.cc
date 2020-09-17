//
// Plugin to readback the TracerProduct.
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

// C++ includes.
#include <iostream>

using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class ReadTracerProduct00 : public art::EDAnalyzer {
  public:
    explicit ReadTracerProduct00(fhicl::ParameterSet const& );
    virtual ~ReadTracerProduct00() { }

    void analyze( art::Event const& e);

  private:

  };

  ReadTracerProduct00::ReadTracerProduct00(fhicl::ParameterSet const& pset )
    : art::EDAnalyzer(pset){
  }

  void
  ReadTracerProduct00::analyze(art::Event const& event) {

    art::Handle<TracerProduct> tpHandle;
    event.getByLabel("tracerTest",tpHandle);

    if ( !tpHandle.isValid() ){
      cout << "Cannot find TracerProduct in the event." << endl;
    }
    TracerProduct const& prod = *tpHandle;

    mf::LogVerbatim("Tracing") << "ReadTracerProduct00:analyze: " << prod;

  } // end of ::analyze.

}

using mu2e::ReadTracerProduct00;
DEFINE_ART_MODULE(ReadTracerProduct00)
