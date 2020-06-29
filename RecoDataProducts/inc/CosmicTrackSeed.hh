#ifndef RecoDataProducts_CosmicTrackSeed_hh
#define RecoDataProducts_CosmicTrackSeed_hh

// Mu2e includes
#include "BTrk/TrkBase/TrkT0.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include <vector>

namespace mu2e {
  
  struct TimeCluster;

  struct CosmicTrackSeed {

    TrkT0 const& t0() const { return _t0; }//from the Time cluster
    CosmicTrack const& track() const { return _track; }
    TrkFitFlag const& status() const { return _status; }
    art::Ptr<TimeCluster> const& timeCluster() const { return _timeCluster; }


    TrkT0	             _t0;	      // t0 for this track
    CosmicTrack              _track;	     // Cosmic track created from these hits  
    TrkFitFlag	             _status;      // status of processes used to create this seed
    // helixOK: ???, helixConverged: ???, circleInit: ???, Straight: ???, hitsOK: ???
    art::Ptr<TimeCluster>    _timeCluster; // associated time cluster
    
    ComboHitCollection       _straw_chits;    // get the straw level hits and store here (need to find panel hits first)
 
    // For future use only
    std::vector<TrkStrawHitSeed> const& trkstrawhits() const { return _trkstrawhits;}
    std::vector<TrkStrawHitSeed>  _trkstrawhits; //vector of associated trkstrawhits
  };
   typedef std::vector<mu2e::CosmicTrackSeed> CosmicTrackSeedCollection;
} 

#endif /* RecoDataProducts_CosmicTrackSeed_hh */
