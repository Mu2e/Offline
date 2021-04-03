//
// A PhysicsConstructor that adds step limiters to some particles.
//
//
// Original author Rob Kutschke
//
// Notes
// 1) These methods are simple enough that one might be tempted to put
//    them in the .hh file.  But that will not work since they are never
//    called by code that we compile; except for the constructor, are
//    all methods callbacks that are called be G4.

// Mu2e includes
#include "Mu2eG4/inc/Mu2eG4StepLimiterPhysicsConstructor.hh"
#include "Mu2eG4/inc/addStepLimiter.hh"

using namespace std;

namespace mu2e {

  Mu2eG4StepLimiterPhysicsConstructor::Mu2eG4StepLimiterPhysicsConstructor():
    G4VPhysicsConstructor("Mu2eG4StepLimiterPhysicsConstructor"){
  }

  Mu2eG4StepLimiterPhysicsConstructor::~Mu2eG4StepLimiterPhysicsConstructor(){
  }

  void Mu2eG4StepLimiterPhysicsConstructor::ConstructParticle(){
  }

  void Mu2eG4StepLimiterPhysicsConstructor::ConstructProcess(){
    addStepLimiter();
  }


}  // end namespace mu2e
