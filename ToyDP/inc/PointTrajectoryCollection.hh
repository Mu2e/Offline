#ifndef ToyDP_PointTrajectoryCollection_hh
#define ToyDP_PointTrajectoryCollection_hh

//
// Define a type for a collection of PointTrajectory objects.
// The key is the simulated particle ID (same as for the
// SimParticleCollection).
//
// $Id: PointTrajectoryCollection.hh,v 1.3 2011/05/18 02:27:19 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:19 $
//
// Original author Rob Kutschke
//

#include "GeneralUtilities/inc/MapVector.hh"
#include "ToyDP/inc/PointTrajectory.hh"

namespace mu2e {
   typedef MapVector<mu2e::PointTrajectory> PointTrajectoryCollection;
}

#endif /* ToyDP_PointTrajectoryCollection_hh */
