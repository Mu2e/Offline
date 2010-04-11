#ifndef AddStepLimiter_HH
#define AddStepLimiter_HH
//
// Free functions to add step limiters to some particles.
//
// $Id: addStepLimiter.hh,v 1.1 2010/04/11 15:15:12 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/11 15:15:12 $
//
// Original author Rob Kutschke
//
// Notes:
// 1) Getting G4 to limit step sizes requires two steps.
//     a) Create a step limiter on a per particle-type basiss.
//     b) Add a user step limit for those logical volumes
//        for which the step limiter is to be active.
//    This routine only does part a).

#include <vector>

#include "G4String.hh"

namespace mu2e{

  // Add step limiters for a standard list of particles:
  // e, mu, pi, K, p (particles and anti-particles) plus chargedgeantino.
  void addStepLimiter ();

  // Add step limiters for a user specified list of particles.
  void addStepLimiter ( const std::vector<G4String>& list);

}  // end namespace mu2e

#endif


