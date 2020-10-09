#ifndef Mu2eG4_physicsListDecider_hh
#define Mu2eG4_physicsListDecider_hh
//
// Decide which physics list to use.
//
//
// Original author Rob Kutschke
//

#include "Mu2eG4/inc/Mu2eG4Config.hh"

// Forward declarations
class G4VUserPhysicsList;

namespace mu2e{

  // The returned pointer should be passed to G4, which will
  // take ownership of it.
  G4VUserPhysicsList* physicsListDecider (const Mu2eG4Config::Physics& phys, const Mu2eG4Config::Debug& debug);


}  // end namespace mu2e
#endif /* Mu2eG4_physicsListDecider_hh */
