#ifndef Mu2eG4_addPointTrajectories_hh
#define Mu2eG4_addPointTrajectories_hh
//
// Dig the trajectories out of the G4 internals and add them to the event.
//
// $Id: addPointTrajectories.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
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

#endif /* Mu2eG4_addPointTrajectories_hh */


