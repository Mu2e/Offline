// Andrei Gaponenko, 2014

#ifndef MCDataProducts_inc_SimParticleTrackerHitSummaryAssns_hh
#define MCDataProducts_inc_SimParticleTrackerHitSummaryAssns_hh

#include "art/Persistency/Common/Assns.h"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleTrackerHitSummaryCollection.hh"

namespace mu2e {
  // Used as a one-to-one assn
  typedef art::Assns<SimParticle,SimParticleTrackerHitSummary> SimParticleTrackerHitSummaryAssns;
}

#endif /* MCDataProducts_inc_SimParticleTrackerHitSummaryAssns_hh */
