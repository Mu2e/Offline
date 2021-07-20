#ifndef Mu2eG4_Mu2eG4MinimalModularPhysicsList_hh
#define Mu2eG4_Mu2eG4MinimalModularPhysicsList_hh
//
// Define a minimal physics list for G4, just transportation.
// Used for debugging geometry.
//
//
// Original author Rob Kutschke
//

#include "Geant4/G4VModularPhysicsList.hh"

namespace mu2e {
  class Mu2eG4MinimalModularPhysicsList: public G4VModularPhysicsList{
  public:
    Mu2eG4MinimalModularPhysicsList();
    ~Mu2eG4MinimalModularPhysicsList();

  protected:

    // These methods are called by G4 not by users.
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();

  };

}  // end namespace mu2e

#endif /* Mu2eG4_Mu2eG4MinimalModularPhysicsList_hh */


