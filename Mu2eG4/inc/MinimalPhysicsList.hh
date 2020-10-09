#ifndef Mu2eG4_MinimalPhysicsList_hh
#define Mu2eG4_MinimalPhysicsList_hh
//
// Define a minimal physics list for G4, just transportation.
// Used for debugging geometry.
//
//
// Original author Rob Kutschke
//

#include "G4VModularPhysicsList.hh"

namespace mu2e {
  class MinimalPhysicsList: public G4VModularPhysicsList{
  public:
    MinimalPhysicsList();
    ~MinimalPhysicsList();

  protected:

    // These methods are called by G4 not by users.
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();

  };

}  // end namespace mu2e

#endif /* Mu2eG4_MinimalPhysicsList_hh */


