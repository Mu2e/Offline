//
// Tests for the bad Ptr bug
//
// $Id: PtrBug01_module.cc,v 1.1 2011/11/03 02:01:44 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/11/03 02:01:44 $
//
// Original author Rob Kutschke
//

// Framework includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes.
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"

// C++ includes.
#include <iostream>
#include <string>
#include <set>

using namespace std;

namespace mu2e {

  class PtrBug01 : public art::EDAnalyzer {
  public:

    explicit PtrBug01(fhicl::ParameterSet const& pset);
    virtual ~PtrBug01() { }

    // The framework calls this for each event.
    void analyze(const art::Event& e);

  private:

  };

  PtrBug01::PtrBug01(fhicl::ParameterSet const& pset) {
  }

  void PtrBug01::analyze(const art::Event& event) {

    art::Handle<SimParticleCollection> simsHandle;
    event.getByLabel("g4run", simsHandle);
    SimParticleCollection const& sims(*simsHandle);

    cout << "Sims.size: " << sims.size() << endl;
    for ( SimParticleCollection::const_iterator i=sims.begin(), e=sims.end();
          i != e; ++i){
      SimParticle const& sim(i->second);
      cout << "Sim: "
           << sim.id() << " "
           << sim.genParticle().id()  <<  " "
           << sim.genParticle().key() <<  " "
           << sim.parent().id()       <<  " "
           << sim.parent().key()
           << endl;
      if ( sim.hasParent() ){
        cout << ":   "
             << sim.parent()->id()
             << endl;
      }
    }


  }

}  // end namespace mu2e

// Part of the magic that makes this class a module.
// create an instance of the module.  It also registers
using mu2e::PtrBug01;
DEFINE_ART_MODULE(PtrBug01);
