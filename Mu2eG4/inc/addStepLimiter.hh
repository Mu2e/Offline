#ifndef Mu2eG4_addStepLimiter_hh
#define Mu2eG4_addStepLimiter_hh
//
// Free functions to add step limiters to some particles.
//
// $Id: addStepLimiter.hh,v 1.2 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:36 $
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

#endif /* Mu2eG4_addStepLimiter_hh */


