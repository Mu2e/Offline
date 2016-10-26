// A persistable class to correlate tracks (TrackSummary) with MC particles.
//
// Andrei Gaponenko, 2014

#ifndef MCDataProducts_inc_TrackSummaryTruthAssns_hh
#define MCDataProducts_inc_TrackSummaryTruthAssns_hh

#include <ostream>

#include "canvas/Persistency/Common/Assns.h"

#include "MCDataProducts/inc/SimParticle.hh"
#include "RecoDataProducts/inc/TrackSummary.hh"

namespace mu2e {

  //================================================================
  // Information attached to track-truth pair: counts hits
  // attached to both the track and the particle.

  class TrackSummaryMatchInfo {
  public:
    TrackSummaryMatchInfo(unsigned nPrincipal, unsigned nAll)
      : nPrincipal_(nPrincipal), nAll_(nAll)
    {}

    // Number of active hits on track that were "made" by the particle
    unsigned nPrincipal()   const { return nPrincipal_; }

    // Number of active hits on track to which the particle contributed
    unsigned nAll()    const { return nAll_; }

    // for persistency
    TrackSummaryMatchInfo() : nPrincipal_(-1u), nAll_(-1u) {}

  private:
    unsigned nPrincipal_;
    unsigned nAll_;
  };

  // This is a many-to-many Assns, in general
  typedef art::Assns<SimParticle,TrackSummary,TrackSummaryMatchInfo> TrackSummaryTruthAssns;

  std::ostream& operator<<(std::ostream& os, const TrackSummaryMatchInfo& mi);

} // namespace mu2e

#endif /* MCDataProducts_inc_TrackSummaryTruthAssns_hh */
