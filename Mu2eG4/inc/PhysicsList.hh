#ifndef PhysicsList_HH
#define PhysicsList_HH
//
// Construct particles; construct and register physics processes.
//
// This is modeled on:
//   $G4INSTALL/examples/novice/N02/include/ExN02PhysicsList.hh 
//    with cvs tag: version 1.12 2008/09/22 16:41:20 maire
//
//
// $Id: PhysicsList.hh,v 1.3 2010/04/11 15:16:55 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/11 15:16:55 $
//
// Original author Rob Kutschke
//

// G4 includes
#include "G4VUserPhysicsList.hh"
#include "globals.hh"

namespace mu2e {

  // Forward class references.
  class SimpleConfig;

  class PhysicsList: public G4VUserPhysicsList{

  public:
    PhysicsList( const SimpleConfig& config);
    ~PhysicsList();
    
  protected:

    // Methods called by G4 (should not be called by anyone else).
    void ConstructParticle();
    void ConstructProcess(); 
    void SetCuts();
    
    // Methods used to break up ConstructParticle into smaller pieces.
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBaryons();
    void ConstructAllOthers();
  
    // Methods used to break up ConstructProcess into smaller pieces.
    void ConstructGeneral();
    void ConstructEM();

  private:

    // Run time configuration.
    const SimpleConfig* _config;

  };  

}  // end namespace mu2e

#endif
