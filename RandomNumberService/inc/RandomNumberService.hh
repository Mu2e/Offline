#ifndef RandomNumberService_RandomNumberService_hh
#define RandomNumberService_RandomNumberService_hh
//
// Maintain multiple independent chains of random numbers,
// including save and restore of seeds.  For now it only
// knows about the CLHEP global instances.
//
// $Id: RandomNumberService.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//
//  The present implementation just seeds the CLHEP global
//  random number generator, HepRandom.
//
//  In a later implementation it will:
//   1) Propagate the state of the CLHEP global random number
//      generator across jobs.
//   2) Manage (potentially) many other random number streams that
//      can be used by any module.
//  

// C++ include files
#include <string>
#include <memory>

// Framework include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"

namespace mu2e {

  class RandomNumberService {
  public:
    RandomNumberService(const edm::ParameterSet&, edm::ActivityRegistry&);
    ~RandomNumberService();
    
    void postBeginJob();
    void postEndJob();

  private:
    
    // Seed for the HepRandom singleton.
    long _globalSeed;
    
  };
}

#endif
