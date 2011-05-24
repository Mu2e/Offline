#ifndef Mu2eG4_addPointTrajectories_hh
#define Mu2eG4_addPointTrajectories_hh
//
// Dig the trajectories out of the G4 internals and add them to the event.
//
// $Id: addPointTrajectories.hh,v 1.3 2011/05/24 17:19:03 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:19:03 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

class G4Event;

namespace mu2e{

  // Public entry point.
  void addPointTrajectories ( const G4Event*             event,
                              PointTrajectoryCollection& hits,
                              CLHEP::Hep3Vector const&   mu2eOriginInWorld );

}  // end namespace mu2e

#endif /* Mu2eG4_addPointTrajectories_hh */


