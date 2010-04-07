#ifndef PhysicsListDecider_HH
#define PhysicsListDecider_HH
//
// Decide which physics list to use.
//
// $Id: physicsListDecider.hh,v 1.1 2010/04/07 22:08:58 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/07 22:08:58 $
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
#endif


