#ifndef Mu2eG4_physicsListDecider_hh
#define Mu2eG4_physicsListDecider_hh
//
// Decide which physics list to use.
//
// $Id: physicsListDecider.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
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
