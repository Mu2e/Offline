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
  
  class TimeCluster;

  struct CosmicTrackSeed {

    TrkT0 const& t0() const { return _t0; }//from the Time cluster
    const std::vector<StrawHitIndex>& _shits       () const { return _strawHitIdxs; }
    std::vector<TrkStrawHitSeed> const& trkstrawhits() const { return _trkstrawhits;}
    ComboHitCollection const& hits() const { return _panel_hits; }
    CosmicTrack const& track() const { return _track; }
    TrkFitFlag const& status() const { return _status; }
    art::Ptr<TimeCluster> const& timeCluster() const { return _timeCluster; }
/*///IN DEVELOPMENT: NEW INFRASTRUCTURE////////
    typedef art::Ptr<ComboHit>          ComboHitPtr;
    typedef std::vector<ComboHitPtr>    ComboHitPtrVector;
////////////////////////////*/
    TrkT0	             _t0;	      // t0 for this track
    ComboHitCollection       _panel_hits;	      // hits for track (panel hits)
    ComboHitCollection       _straw_chits;    // get the straw level hits and store here (need to find panel hits first)
   
/*/IN DEVELOPMENT: NEW INFRASTRUCTURE/////
    ComboHitPtrVector	     _panelHits;
    ComboHitPtrVector	     _strawHits;
//////////////////////////*/
    CosmicTrack              _track;	     // Cosmic track created from these hits  
   
    TrkFitFlag	             _status;      // status of processes used to create this seed
    art::Ptr<TimeCluster>    _timeCluster; // associated time cluster
    std::vector<StrawHitIndex> _strawHitIdxs; // associated straw hits: can be empty
    std::vector<TrkStrawHitSeed>  _trkstrawhits; //vector of associated trkstrawhits
   
   
  };
   typedef std::vector<mu2e::CosmicTrackSeed> CosmicTrackSeedCollection;
} 

#endif /* RecoDataProducts_CosmicTrackSeed_hh */
