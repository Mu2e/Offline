#ifndef MCDataProducts_PointTrajectoryCollection_hh
#define MCDataProducts_PointTrajectoryCollection_hh

//
// Define a type for a collection of PointTrajectory objects.
// The key is the simulated particle ID (same as for the
// SimParticleCollection).
//
// $Id: PointTrajectoryCollection.hh,v 1.1 2011/05/24 17:16:43 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/05/24 17:16:43 $
//
// Original author Rob Kutschke
//

#include "GeneralUtilities/inc/MapVector.hh"
#include "MCDataProducts/inc/PointTrajectory.hh"

namespace mu2e {
   typedef MapVector<mu2e::PointTrajectory> PointTrajectoryCollection;
}

#endif /* MCDataProducts_PointTrajectoryCollection_hh */
