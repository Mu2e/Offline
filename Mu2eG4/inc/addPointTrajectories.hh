#ifndef AddPointTrajectories_HH
#define AddPointTrajectories_HH
//
// Dig the trajectories out of the G4 internals and add them to the event.
//
// $Id: addPointTrajectories.hh,v 1.1 2010/11/18 07:22:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/11/18 07:22:43 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "ToyDP/inc/PointTrajectoryCollection.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

class G4Event;

namespace mu2e{

  // Public entry point.
  void addPointTrajectories ( const G4Event*             event,
                              PointTrajectoryCollection& hits,
                              CLHEP::Hep3Vector const&   mu2eOriginInWorld );

}  // end namespace mu2e

#endif


