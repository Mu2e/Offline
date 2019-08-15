#ifndef RecoDataProducts_StraightTrackSeed_hh
#define RecoDataProducts_StraightTrackSeed_hh
//
// S. Middleton - seed class for striaght track fit
// This class is mostly obsolete now-upgrade to cosmictrack

// Mu2e includes
#include "BTrk/TrkBase/TrkT0.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StraightTrack.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include <vector>

namespace mu2e {
  
  class TimeCluster;

  struct StraightTrackSeed {

    ComboHitCollection const& hits() const { return _thits; }
    StraightTrack const& track() const { return _track; }
    TrkFitFlag const& status() const { return _status; }
    art::Ptr<TimeCluster> const& timeCluster() const { return _timeCluster; }

    TrkT0	             _t0;	      // t0 for this track
    ComboHitCollection       _thits;	      // hits for track
    StraightTrack            _track;	     // straight track created from these hits
    TrkFitFlag	             _status;      // status of processes used to create this seed
    art::Ptr<TimeCluster>    _timeCluster; // associated time cluster
  };
   typedef std::vector<mu2e::StraightTrackSeed> StraightTrackSeedCollection;
} // namespace mu2e

#endif /* RecoDataProducts_StraightTrackSeed_hh */
