#ifndef RecoDataProducts_CosmicTrackSeed_hh
#define RecoDataProducts_CosmicTrackSeed_hh
//
// S. Middleton - seed class for striaght track fit
//

// Mu2e includes
#include "BTrk/TrkBase/TrkT0.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include <vector>

namespace mu2e {
  
  class TimeCluster;

  struct CosmicTrackSeed {

    TrkT0 const& t0() const { return _t0; }//from the Time cluster
   // const std::vector<StrawHitIndex>& hits       () const { return _strawHitIdxs; }
    //ComboHitCollection const& hits() const { return _thits; }
    CosmicTrack const& track() const { return _track; }
    TrkFitFlag const& status() const { return _status; }
    art::Ptr<TimeCluster> const& timeCluster() const { return _timeCluster; }

    TrkT0	             _t0;	      // t0 for this track
    ComboHitCollection       _thits;	      // hits for track
    CosmicTrack              _track;	     // Cosmic track created from these hits
    TrkFitFlag	             _status;      // status of processes used to create this seed
    art::Ptr<TimeCluster>    _timeCluster; // associated time cluster
  };
   typedef std::vector<mu2e::CosmicTrackSeed> CosmicTrackSeedCollection;
} // namespace mu2e

#endif /* RecoDataProducts_CosmicTrackSeed_hh */
