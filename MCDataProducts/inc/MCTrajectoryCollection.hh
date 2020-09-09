#ifndef MCDataProducts_MCTrajectoryCollection_hh
#define MCDataProducts_MCTrajectoryCollection_hh

//
// Define a type for a collection of MCTrajectory objects.
//
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
