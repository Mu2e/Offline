//
// Plugin to readback the std::vector<TracerProduct>.
//
// $Id: ReadVectorTracerProduct00_module.cc,v 1.1 2011/06/05 17:30:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/05 17:30:43 $
//
// Original author Rob Kutschke.
//

// Mu2e includes.
#include "Sandbox/inc/TracerProduct.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes.
#include <iostream>
#include <vector>

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

    art::Handle<std::vector<TracerProduct> > tpHandle;
    event.getByLabel("tracerTest",tpHandle);

    if ( !tpHandle.isValid() ){
      cout << "Cannot find std::vector<TracerProduct> in the event." << endl;
    }
    std::vector<TracerProduct> const& prod = *tpHandle;

    for ( size_t i=0; i<prod.size(); ++i){
      mf::LogVerbatim("Tracing") << "      ReadTracerProduct00:analyze: " << prod.at(i);
    }

  } // end of ::analyze.

}

using mu2e::ReadTracerProduct00;
DEFINE_ART_MODULE(ReadTracerProduct00);
