#ifndef Mu2eG4_Mu2eG4MinimalModularPhysicsList_hh
#define Mu2eG4_Mu2eG4MinimalModularPhysicsList_hh
//
// Define a minimal physics list for Geant4
// Used for debugging geometry etc.
//
//
// Original author Rob Kutschke
//

#include "Geant4/G4VModularPhysicsList.hh"

namespace mu2e {
  class Mu2eG4MinimalModularPhysicsList: public G4VModularPhysicsList{
  public:
    Mu2eG4MinimalModularPhysicsList();
    virtual ~Mu2eG4MinimalModularPhysicsList() = default;
    Mu2eG4MinimalModularPhysicsList(const Mu2eG4MinimalModularPhysicsList &) = delete;
    Mu2eG4MinimalModularPhysicsList & operator=(const Mu2eG4MinimalModularPhysicsList &) = delete;
    Mu2eG4MinimalModularPhysicsList & operator=( Mu2eG4MinimalModularPhysicsList && ) = delete;
  };

}  // end namespace mu2e

#endif /* Mu2eG4_Mu2eG4MinimalModularPhysicsList_hh */
