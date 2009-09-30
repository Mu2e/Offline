//
// Maintain multiple independent chains of random numbers,
// including save and restore of seeds.  For now it only
// knows about the CLHEP global instances.
//
// $Id: RandomNumberService.cc,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

// C++ include files
#include <iostream>
#include <typeinfo>

// Framework include files
#include "DataFormats/Provenance/interface/ModuleDescription.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "DataFormats/Provenance/interface/Timestamp.h"
#include "DataFormats/Provenance/interface/LuminosityBlockID.h"
#include "DataFormats/Provenance/interface/RunID.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e include files
#include "RandomNumberService/inc/RandomNumberService.hh"

// Other external includes.
#include "CLHEP/Random/Random.h"

using namespace std;

namespace mu2e {

  RandomNumberService::RandomNumberService(edm::ParameterSet const& iPS, 
					   edm::ActivityRegistry&iRegistry) :
    _globalSeed(iPS.getUntrackedParameter<int>("globalSeed",-1))
  {
    iRegistry.watchPostBeginJob(this, &RandomNumberService::postBeginJob);
    iRegistry.watchPostEndJob(this, &RandomNumberService::postEndJob);
  }
  
  RandomNumberService::~RandomNumberService(){
  }
  
  // This is called after all modules have had their beginJob methods called
  // but before the event loop starts.  I would prefer a preBeginJob
  // method but the corresponding watch method does not exist.
  void 
  RandomNumberService::postBeginJob(){

    // The default random number generator is HepJamesRandom.
    // A comment in the sources says that "seed" should be within 
    // the range [0,900000000]
    // This does not change the full state information - just one seed.
    if ( _globalSeed > -1 ){
      CLHEP::HepRandom::setTheSeed(_globalSeed);
    }
    
  }

  // Later, we can save the state in this method.
  void RandomNumberService::postEndJob(){
  }

  
} // end namespace mu2e
