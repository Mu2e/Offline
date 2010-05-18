//
// Test accessing the RandomNumberService.
//
// $Id: RNGTest_plugin.cc,v 1.4 2010/05/18 21:16:41 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/05/18 21:16:41 $
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
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"

using namespace std;

namespace mu2e {

  class RNGTest : public edm::EDAnalyzer {
  public:
    explicit RNGTest(edm::ParameterSet const& pset):
      _doSkip(pset.getUntrackedParameter<int>("doSkip",0)){
    }
    virtual ~RNGTest() { }

    void analyze( edm::Event const& e, edm::EventSetup const&);

  private:

    int _doSkip;

  };

  void
  RNGTest::analyze( edm::Event const& event, edm::EventSetup const&) {


    // Access Random Number Service.
    edm::Service<edm::RandomNumberGenerator> rng;
    CLHEP::HepRandomEngine& engine = rng->getEngine();

    // Number of calls to this method.
    static int ncalls(-1);

    // Skip some events if requested.
    bool skip = false;
    if ( _doSkip == 1 ){
      if ( (++ncalls)%2 == 0 ) skip = true;
    } else if (_doSkip == 2 ) {
      if ( (++ncalls)%2 == 1 ) skip = true;
    }
    if ( skip ) {
      cerr << "Skipping this event: " << ncalls << endl;
      return;
    }
    
    // Create some distributions that use that engine.
    static CLHEP::RandFlat   flat( engine );
    static CLHEP::RandGaussQ gauss( engine, 10., 2. );

    // Check output of the distributions.
    for ( int i=0; i<5; ++i ){
      cerr << "Fire: "
           << setw(3) << i <<  " "
           << flat.fire()  <<  " "
           << gauss.fire()
           << endl;
    }



  }
  
} // end namespace mu2e

using mu2e::RNGTest;
DEFINE_FWK_MODULE(RNGTest);
