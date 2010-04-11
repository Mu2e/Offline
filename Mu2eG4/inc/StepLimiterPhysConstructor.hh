#ifndef StepLimiterPhysConstructor_HH
#define StepLimiterPhysConstructor_HH
//
// A Physics constructor that adds step limiters to some particles.
//
// $Id: StepLimiterPhysConstructor.hh,v 1.1 2010/04/11 15:15:12 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/11 15:15:12 $
//
// Original author Rob Kutschke
//
#include "G4VPhysicsConstructor.hh"

namespace mu2e {

  class  StepLimiterPhysConstructor: public G4VPhysicsConstructor {

  public:
    StepLimiterPhysConstructor();
    ~StepLimiterPhysConstructor();
    
    void ConstructParticle();
    void ConstructProcess();

  };
  
} // end namespace mu2e

#endif
