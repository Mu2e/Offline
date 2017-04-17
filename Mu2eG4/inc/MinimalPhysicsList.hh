#ifndef Mu2eG4_MinimalPhysicsList_hh
#define Mu2eG4_MinimalPhysicsList_hh
//
// Define a minimal physics list for G4, just transportation.
// Used for debugging geometry.
//
// $Id: MinimalPhysicsList.hh,v 1.3 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
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


