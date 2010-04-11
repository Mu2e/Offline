//
// A PhysicsConstructor that adds step limiters to some particles.
//
// $Id: StepLimiterPhysConstructor.cc,v 1.1 2010/04/11 15:15:12 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/11 15:15:12 $
//
// Original author Rob Kutschke
//
// Notes
// 1) These methods are simple enough that one might be tempted to put
//    them in the .hh file.  But that will not work since they are never
//    called by code that we compile; except for the constructor, are 
//    all methods callbacks that are called be G4.

// Mu2e includes
#include "Mu2eG4/inc/StepLimiterPhysConstructor.hh"
#include "Mu2eG4/inc/addStepLimiter.hh"

using namespace std;

namespace mu2e {

  StepLimiterPhysConstructor::StepLimiterPhysConstructor():
    G4VPhysicsConstructor("StepLimiterPhysConstructor"){
  }

  StepLimiterPhysConstructor::~StepLimiterPhysConstructor(){
  }

  void StepLimiterPhysConstructor::ConstructParticle(){
  }

  void StepLimiterPhysConstructor::ConstructProcess(){
    addStepLimiter();
  }


}  // end namespace mu2e
