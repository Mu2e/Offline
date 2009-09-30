#ifndef PhysicsList_h
#define PhysicsList_h 1
//
// Define the physics list with G4 for the Mu2e.
//
// $Id: PhysicsList.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//

#include "G4VUserPhysicsList.hh"

namespace mu2e {
  class PhysicsList: public G4VUserPhysicsList{
  public:
    PhysicsList();
    ~PhysicsList();

  protected:

    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();

    void SetCuts();
  };

}  // end namespace mu2e

#endif

 
