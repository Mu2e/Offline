//
// Plugin to readback the TracerProductCollection
//
//
// Original author Rob Kutschke.
//

// Mu2e includes.
#include "Sandbox/inc/TracerProduct.hh"
#include "Sandbox/inc/TracerProductCollection.hh"

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

// Other includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// C++ includes.

using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class ReadVectorTracerProduct00 : public art::EDAnalyzer {
  public:
    explicit ReadVectorTracerProduct00(fhicl::ParameterSet const& pset );
    virtual ~ReadVectorTracerProduct00() { }

    void analyze( art::Event const& e);

  private:

  };

  ReadVectorTracerProduct00::ReadVectorTracerProduct00(fhicl::ParameterSet const& pset)
    : art::EDAnalyzer(pset){
  }

  void
  ReadVectorTracerProduct00::analyze(art::Event const& event) {

    art::Handle<TracerProductCollection> tpHandle;
    event.getByLabel("tracerTest",tpHandle);
    TracerProductCollection const& prod = *tpHandle;

    for ( size_t i=0; i<prod.size(); ++i){
      TracerProduct const& p = prod.at(i);
      mf::LogVerbatim("Tracing") << "      ReadVectorTracerProduct00:analyze: " << p << endl;
    }

  } // end of ::analyze.

}

using mu2e::ReadVectorTracerProduct00;
DEFINE_ART_MODULE(ReadVectorTracerProduct00)
