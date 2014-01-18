#ifndef MCDataProducts_MCTrajectoryCollection_hh
#define MCDataProducts_MCTrajectoryCollection_hh

//
// Define a type for a collection of MCTrajectory objects.
//
// $Id: MCTrajectoryCollection.hh,v 1.1 2014/01/18 03:08:28 kutschke Exp $
// $Author: kutschke $
// $Date: 2014/01/18 03:08:28 $
//
// Contact person Rob Kutschke
//

#include "MCDataProducts/inc/MCTrajectory.hh"
#include "MCDataProducts/inc/SimParticle.hh"

#include <map>

namespace mu2e {
  typedef std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory> MCTrajectoryCollection;
}

#endif /* MCDataProducts_MCTrajectoryCollection_hh */
