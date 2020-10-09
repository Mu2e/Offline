//
// Plugin to readback the std::vector<TracerProduct>.
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
#include <vector>

using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class ReadVectorTracerProduct00 : public art::EDAnalyzer {
  public:
    explicit ReadVectorTracerProduct00(fhicl::ParameterSet const& );
    virtual ~ReadVectorTracerProduct00() { }

    void analyze( art::Event const& e);

  private:

  };

  ReadVectorTracerProduct00::ReadVectorTracerProduct00(fhicl::ParameterSet const& pset)
    : art::EDAnalyzer(pset){
  }


  void
  ReadVectorTracerProduct00::analyze(art::Event const& event) {

    art::Handle<std::vector<TracerProduct> > tpHandle;
    event.getByLabel("tracerTest",tpHandle);

    if ( !tpHandle.isValid() ){
      cout << "Cannot find std::vector<TracerProduct> in the event." << endl;
    }
    std::vector<TracerProduct> const& prod = *tpHandle;

    for ( size_t i=0; i<prod.size(); ++i){
      mf::LogVerbatim("Tracing") << "      ReadVectorTracerProduct00:analyze: " << prod.at(i);
    }

  } // end of ::analyze.

}

using mu2e::ReadVectorTracerProduct00;
DEFINE_ART_MODULE(ReadVectorTracerProduct00)
