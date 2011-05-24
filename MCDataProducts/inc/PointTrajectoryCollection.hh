#ifndef MCDataProducts_PointTrajectoryCollection_hh
#define MCDataProducts_PointTrajectoryCollection_hh

//
// Define a type for a collection of PointTrajectory objects.
// The key is the simulated particle ID (same as for the
// SimParticleCollection).
//
// $Id: PointTrajectoryCollection.hh,v 1.2 2011/05/24 20:03:31 wb Exp $
// $Author: wb $
// $Date: 2011/05/24 20:03:31 $
//
// Original author Rob Kutschke
//

#include "MCDataProducts/inc/PointTrajectory.hh"
#include "cetlib/map_vector.h"

namespace mu2e {
   typedef cet::map_vector<mu2e::PointTrajectory> PointTrajectoryCollection;
}

#endif /* MCDataProducts_PointTrajectoryCollection_hh */
