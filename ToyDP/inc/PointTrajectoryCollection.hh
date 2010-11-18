#ifndef ToyDP_PointTrajectoryCollection_hh
#define ToyDP_PointTrajectoryCollection_hh

//
// Define a type for a collection of PointTrajectory objects.
// The key is the simulated particle ID (same as for the 
// SimParticleCollection).
//
// $Id: PointTrajectoryCollection.hh,v 1.1 2010/11/18 07:21:39 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/11/18 07:21:39 $
//
// Original author Rob Kutschke
//

#include "GeneralUtilities/inc/MapVector.hh"
#include "ToyDP/inc/PointTrajectory.hh"

namespace mu2e {
   typedef MapVector<mu2e::PointTrajectory> PointTrajectoryCollection;
}

#endif
