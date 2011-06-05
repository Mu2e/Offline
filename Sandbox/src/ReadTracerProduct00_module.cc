//
// Plugin to readback the trancer data product.
//
// $Id: ReadTracerProduct00_module.cc,v 1.1 2011/06/05 16:41:14 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/05 16:41:14 $
//
// Original author Rob Kutschke.
//

// Mu2e includes.
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "Sandbox/inc/TracerProduct.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
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
    explicit ReadTracerProduct00(fhicl::ParameterSet const& pset){}
    virtual ~ReadTracerProduct00() { }

    void analyze( art::Event const& e);

  private:

  };

  void
  ReadTracerProduct00::analyze(art::Event const& event) {

    art::Handle<TracerProduct> tpHandle;
    event.getByLabel("tracerTest",tpHandle);

    if ( !tpHandle.isValid() ){
      cout << "Cannot find TracerProduct00Collection in the event." << endl;
    }
    TracerProduct const& prod = *tpHandle;

    mf::LogVerbatim("Tracing") << "ReadTracerProduct00:analyze: " << prod;

  } // end of ::analyze.

}

using mu2e::ReadTracerProduct00;
DEFINE_ART_MODULE(ReadTracerProduct00);
