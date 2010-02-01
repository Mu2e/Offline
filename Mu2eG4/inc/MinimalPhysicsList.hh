#ifndef MinimalPhysicsList_h
#define MinimalPhysicsList_h 1
//
// Define a minimal physics list for G4, just transportation.
// Used for debugging geometry.
//
// $Id: MinimalPhysicsList.hh,v 1.1 2010/02/01 00:16:39 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/02/01 00:16:39 $
//
// Original author Rob Kutschke
//

#include "G4VUserPhysicsList.hh"

namespace mu2e {
  class MinimalPhysicsList: public G4VUserPhysicsList{
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

#endif

 
