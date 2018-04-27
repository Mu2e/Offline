//
// Tests for the bad Ptr bug
//
// $Id: PtrBug01_module.cc,v 1.5 2013/10/21 21:01:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/10/21 21:01:23 $
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
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"

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

    int which_;

    void bug01a( art::Event const& event);
    void bug01b( art::Event const& event);
    void bug01c( art::Event const& event);
    void bug01d( art::Event const& event);
    void bug01e( art::Event const& event);
    void bug01f( art::Event const& event);


  };

  PtrBug01::PtrBug01(fhicl::ParameterSet const& pset )
    : art::EDAnalyzer(pset){
    cerr << "Select a test by entering a number: " << endl;
    cerr << "   1 - From std::vector into another std::vector." << endl;
    cerr << "   2 - From std::vector into a map_vector." << endl;
    cerr << "   3 - From mapvector into a std::vector" << endl;
    cerr << "   4 - Self reference into the same mapvector, mother" << endl;
    cerr << "   5 - Self reference into the same mapvector, daughters" << endl;
    cerr << "   6 - From vector<ptr> to vector."  << endl;
    cin >> which_;
  }

  void PtrBug01::analyze(const art::Event& event) {

    switch (which_){

    case 1:
      bug01a(event);
      break;

    case 2:
      bug01b(event);
      break;

    case 3:
      bug01c(event);
      break;

    case 4:
      bug01d(event);
      break;

    case 5:
      bug01e(event);
      break;

    case 6:
      bug01f(event);
      break;

    default:
      cerr << "No such test: " << endl;

    }
  }

  // For an object inside a std::vector, follow a Ptr to
  // another object in a different std::vector.
  void PtrBug01::bug01a(const art::Event& event) {

    art::Handle<CaloCrystalHitCollection> crystalsHandle;
    event.getByLabel("CaloCrystalHitsMaker", crystalsHandle);
    CaloCrystalHitCollection const& crystals(*crystalsHandle);
    cout << "CrystalHits.size: " << event.id() << " " <<  crystals.size() << endl;

    typedef CaloCrystalHitCollection::const_iterator Iter;
    for ( Iter i=crystals.begin(), e=crystals.end(); i != e; ++i){
      CaloCrystalHit const& crystal(*i);
      cerr << "Crystal: "
           << crystal.id() << " "
           << crystal.energyDep()  <<  " "
           << crystal.nROId()
           << endl;
      /*
      std::vector<art::Ptr<CaloHit> > const & readouts = crystal.readouts();
      for ( std::vector<art::Ptr<CaloHit> >::const_iterator j=readouts.begin(), je=readouts.end();
            j != je; ++j ){
        art::Ptr<CaloHit> const& p(*j);
        cerr << "   "
             << p.id() << " "
             << p.key() << " "
             << p.productGetter()
             << endl;
        cerr << "      "
             << p->id() << " "
             << p->energyDep() << " "
             << endl;
      }
      */
    }

  } // end bug01a

  // For an object inside a std::vector, follow a Ptr to
  // an object in a map_vector.
  void PtrBug01::bug01b(const art::Event& event) {

    art::Handle<StepPointMCCollection> stepsHandle;
    event.getByLabel("g4run", "tracker", stepsHandle);
    StepPointMCCollection const& steps(*stepsHandle);

    cout << "Steps.size: " << event.id() << " " << steps.size() << endl;
    for ( StepPointMCCollection::const_iterator i=steps.begin(), e=steps.end();
          i != e; ++i){
      StepPointMC const& step(*i);

      cout << "Step: "
           << step.volumeId()          << " "
           << step.totalEDep()*1000.   <<  " "  // keV
           << step.simParticle().id()  <<  " "
           << step.simParticle().key() <<  " | "
           << step.simParticle().productGetter()
           << endl;

      cout << "   Step 1: "
           << step.simParticle()->id()
           << endl;
      cout << "   Step 2: " << endl;
    }

  } // end bug01b

  // For an object inside a map_vector, follow a Ptr to
  // an object in a std::vector.
  void PtrBug01::bug01c(const art::Event& event) {

    art::Handle<SimParticleCollection> simsHandle;
    event.getByLabel("g4run", simsHandle);
    SimParticleCollection const& sims(*simsHandle);

    cout << "Sims.size: " << event.id() << " " << sims.size() << endl;
    for ( SimParticleCollection::const_iterator i=sims.begin(), e=sims.end();
          i != e; ++i){
      SimParticle const& sim(i->second);
      cout << "Sim: "
           << sim.id() << " "
           << sim.genParticle().id()  <<  " "
           << sim.genParticle().key() <<  " "
           << sim.genParticle().productGetter()
           << endl;
      if ( sim.isPrimary() ){
        cout << "   Gen 1: "
             << sim.genParticle()->pdgId()
             << endl;
        cout << "   Gen 2: " << endl;
      }
    }

  } // end bug01c

  // For an object inside a map_vector, follow a Ptr to
  // another object in the same map_vector.
  void PtrBug01::bug01d(const art::Event& event) {

    art::Handle<SimParticleCollection> simsHandle;
    event.getByLabel("g4run", simsHandle);
    SimParticleCollection const& sims(*simsHandle);

    cout << "Sims.size: " << sims.size() << endl;
    for ( SimParticleCollection::const_iterator i=sims.begin(), e=sims.end();
          i != e; ++i){
      SimParticle const& sim(i->second);
      cout << "Sim: "
           << sim.id()           << " "
           << sim.parent().id()  <<  " "
           << sim.parent().key() <<  " "
           << sim.parent().productGetter()
           << endl;
      if ( sim.hasParent() ){
        cout << ":   "
             << sim.parent()->id()
             << endl;
      }
    }

  } // end bug01d

  // For an object inside a map_vector, follow a vector<Ptr> to
  // objects in the same map_vector.
  void PtrBug01::bug01e(const art::Event& event) {

    art::Handle<SimParticleCollection> simsHandle;
    event.getByLabel("g4run", simsHandle);
    SimParticleCollection const& sims(*simsHandle);

    cout << "Sims.size: " << sims.size() << endl;
    for ( SimParticleCollection::const_iterator i=sims.begin(), e=sims.end();
          i != e; ++i){
      SimParticle const& sim(i->second);
      std::vector<art::Ptr<SimParticle> > const& daughters(sim.daughters());
      cout << "Sim: "
           << sim.id()           << " "
           << daughters.size()   <<  " "
           << endl;
      for ( std::vector<art::Ptr<SimParticle> >::const_iterator j=daughters.begin(), je=daughters.end();
            j != je; ++j){
        art::Ptr<SimParticle> const& dau(*j);
        cout << "  Dau: "
             << dau.id() << " "
             << dau.key() << " "
             << dau.productGetter() << " "
             << endl;
        cout << "          " <<  dau->id() << endl;
      }
    }

  } // end bug01e

  // For an object inside a map_vector, follow a vector<Ptr> to
  // objects in the same map_vector.
  void PtrBug01::bug01f(const art::Event& event) {

    art::Handle<StrawHitCollection> strawsHandle;
    event.getByLabel("makeSH", strawsHandle);
    StrawHitCollection const& straws(*strawsHandle);

    art::Handle<PtrStepPointMCVectorCollection> stepPtrsHandle;
    event.getByLabel("makeSH", "StrawHitMCPtr", stepPtrsHandle);
    PtrStepPointMCVectorCollection const& stepPtrs(*stepPtrsHandle);

    cout << "Straws.size: " << straws.size() << endl;
    for ( size_t i=0; i<straws.size(); ++i){
      StrawHit const& straw(straws[i]);
      PtrStepPointMCVector const& steps(stepPtrs[i]);
      cout << "Straw: "
           << straw.strawId() << " "
           << straw.energyDep()*1000.<< " " // keV
           << steps.size()
           << endl;
      double sum(0.);
      for ( size_t j=0; j<steps.size(); ++j){
        art::Ptr<StepPointMC> p(steps[j]);
        cout << "    "
             << j << " "
             << p.id() << " "
             << p.key() << " "
             << p.productGetter() << " "
             << endl;
        cout << "          "
             << p->trackId()  << " "
             << p->volumeId() << " "
             << p->totalEDep()*1000. // keV
             << endl;
        sum += p->totalEDep();
      }
      cout << "    Sum: "
           << sum*1000.              << " "
           << straw.energyDep()*1000.<< " " // keV
           << ( sum - straw.energyDep())*1000
           << endl;
    }

  } // end bug01f


}  // end namespace mu2e

// Part of the magic that makes this class a module.
// create an instance of the module.  It also registers
using mu2e::PtrBug01;
DEFINE_ART_MODULE(PtrBug01)
