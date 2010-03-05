//
// Test accessing the RandomNumberService.
//
// $Id: RNGTest_plugin.cc,v 1.1 2010/03/05 16:07:38 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/03/05 16:07:38 $
//
// Original author Rob Kutschke
//

// C++ includes.
#include <iostream>
#include <string>

// Framework includes.
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e includes
#include "RandomNumberService/inc/RandomNumberService.hh"

// CLHEP includes
#include "CLHEP/Random/RandomEngine.h"

using namespace std;

namespace mu2e {

  class RNGTest : public edm::EDAnalyzer {
  public:
    explicit RNGTest(edm::ParameterSet const& pset){}
    virtual ~RNGTest() { }

    void analyze( edm::Event const& e, edm::EventSetup const&);

  private:

  };

  void
  RNGTest::analyze( edm::Event const& event, edm::EventSetup const&) {

    cerr << "\nHave event: " << event.id() << endl;

    edm::Service<RandomNumberService> rns;
    cerr << "Is concrete available: "  << rns.isAvailable() << endl;

    edm::Service<edm::RandomNumberGenerator> rng;
    cerr << "Is base available:     "  << rng.isAvailable() << endl;

    CLHEP::HepRandomEngine& engine = rng->getEngine();
    cerr << "Fire: " << engine.flat() << endl;

    vector<unsigned long> v = engine.put();
    cerr << "Size: " << v.size() << endl;

    rng->print();

    //if ( event.id().event() == 1 ){
    //engine.put(cout);
    //}

  }
  
}


using mu2e::RNGTest;
DEFINE_FWK_MODULE(RNGTest);
