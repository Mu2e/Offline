//
// Plugin to readback the transient data product.
//
//
// Original author Rob Kutschke.
//

// C++ includes.
#include <iostream>

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"

// Mu2e includes.
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "Sandbox/inc/TransientProduct00Collection.hh"

using namespace std;

namespace mu2e {

  //--------------------------------------------------------------------
  //
  //
  class ReadTransientProduct00 : public art::EDAnalyzer {
  public:
    explicit ReadTransientProduct00(fhicl::ParameterSet const& );
    virtual ~ReadTransientProduct00() { }

    void analyze( art::Event const& e);

  private:

  };

  ReadTransientProduct00::ReadTransientProduct00(fhicl::ParameterSet const& pset)
    : art::EDAnalyzer(pset){
  }


  void
  ReadTransientProduct00::analyze(art::Event const& event) {

    art::Handle<StrawHitCollection> hitsHandle;
    event.getByLabel("makeSH",hitsHandle);
    StrawHitCollection const& hits = *hitsHandle;

    art::Handle<TransientProduct00Collection> tpHandle;
    event.getByLabel("transientTest",tpHandle);

    if ( !tpHandle.isValid() ){
      cout << "Cannot find TransientProduct00Collection in the event." << endl;
    }
    TransientProduct00Collection const& prod = *tpHandle;


    if ( prod.size() != hits.size() ) {
      cout << "Error: Read sizes are: "
           << prod.size()  << " "
           << hits.size()  << " "
           << endl;
    }

    if ( hits.size() != prod.size() ) return;

    int nbad(0);
    for ( size_t i=0; i<prod.size(); ++i ){
      StrawHit const& s0 = hits.at(i);
      StrawHit const& s1 = prod.at(i).hit();
      //if ( s0.strawId() != s1.strawId() ) {
      if ( &s0 != &s1 ) {
        ++nbad;
        cout << "Event: "
             << event.id()
             << " Straw: "
             << s0.strawId() << " "
             << s1.strawId() << " "
             << endl;
      }

    }

    if ( nbad > 0 ){
      cout << "Error: Number of bad comparisons: " << nbad << endl;
    } else{
      cout << "All OK." << prod.size() << " " << hits.size() << endl;
    }


  } // end of ::analyze.

}


using mu2e::ReadTransientProduct00;
DEFINE_ART_MODULE(ReadTransientProduct00)
