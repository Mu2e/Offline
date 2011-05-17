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

// Forward declarations
class G4VUserPhysicsList;

namespace mu2e{

  // Forward declarations within Mu2e namespace.
  class SimpleConfig;

  // The returned pointer should be passed to G4, which will
  // take ownership of it.
  G4VUserPhysicsList* physicsListDecider ( const SimpleConfig& config );


}  // end namespace mu2e
#endif /* Mu2eG4_physicsListDecider_hh */


