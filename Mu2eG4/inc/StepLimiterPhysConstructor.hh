#ifndef Mu2eG4_StepLimiterPhysConstructor_hh
#define Mu2eG4_StepLimiterPhysConstructor_hh
//
// A Physics constructor that adds step limiters to some particles.
//
// $Id: StepLimiterPhysConstructor.hh,v 1.3 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
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

#endif /* Mu2eG4_StepLimiterPhysConstructor_hh */
