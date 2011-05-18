#ifndef Mu2eG4_PhysicsList_hh
#define Mu2eG4_PhysicsList_hh
//
// Construct particles; construct and register physics processes.
//
// This is modeled on:
//   $G4INSTALL/examples/novice/N02/include/ExN02PhysicsList.hh
//    with cvs tag: version 1.12 2008/09/22 16:41:20 maire
//
//
// $Id: PhysicsList.hh,v 1.5 2011/05/18 02:27:17 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:17 $
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

#endif /* Mu2eG4_PhysicsList_hh */
