#ifndef RecoDataProducts_TrackSeed_hh
#define RecoDataProducts_TrackSeed_hh
//
// out data of the time pattern recognition algorithm to seed the Kalman filter (used for the ITracker)
//
// $Id: $
// $Author:  $
// $Date: $
//
// Original author G. Tassielli
//

// Mu2e includes
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/HelixVal.hh"

namespace mu2e {

  struct TrackSeed{
  
    TimeCluster        _timeCluster;
    HelixVal           _helix;
  
  public:

    TrackSeed() {}
    
    double d0      () const {return _helix.d0();    }
    double phi0    () const {return _helix.phi0();  }
    double omega   () const {return _helix.omega(); }
    double z0      () const {return _helix.z0();    }
    double tanDip  () const {return _helix.tanDip();}

    TrkT0 const& t0      () const {return _timeCluster.t0();    }
    
  };



} // namespace mu2e

#endif /* RecoDataProducts_TrackSeed_hh */
