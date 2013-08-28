#ifndef Mu2eG4_addPointTrajectories_hh
#define Mu2eG4_addPointTrajectories_hh
//
// Dig the trajectories out of the G4 internals and add them to the event.
//
// $Id: addPointTrajectories.hh,v 1.4 2013/08/28 05:58:17 gandr Exp $
// $Author: gandr $
// $Date: 2013/08/28 05:58:17 $
//
// Original author Rob Kutschke
//

// Mu2e includes
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

class G4Event;

namespace mu2e{

  class SimParticleHelper;

  // Public entry point.
  void addPointTrajectories ( const G4Event*             event,
                              PointTrajectoryCollection& hits,
                              const SimParticleHelper& spHelper,
                              CLHEP::Hep3Vector const&   mu2eOriginInWorld );

}  // end namespace mu2e

#endif /* Mu2eG4_addPointTrajectories_hh */


