//
// Functor to calculate a track t0 from a TimeCluster.  This takes into account
// the specific offsets of each subsystem
//
// $Id: HelixFit.hh,v 1.8 2014/07/10 14:47:26 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/07/10 14:47:26 $
//
#ifndef TrkReco_TrkTimeCalculator_HH
#define TrkReco_TrkTimeCalculator_HH

// framework
#include "fhiclcpp/ParameterSet.h"
// data
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/TrkFitDirection.hh"
// HelixFit objects
#include "RecoDataProducts/inc/HelixSeed.hh"
// BaBar
#include "BTrk/TrkBase/TrkErrCode.hh"
//root
class TH1F;
// C+

namespace mu2e 
{
  class TrkTimeCalculator {
    public:
      // parameter set must be passed in on construction
      explicit TrkTimeCalculator(fhicl::ParameterSet const&);
      virtual ~TrkTimeCalculator();
      // update the t0 value inside a TimeCluster
      void updateT0(TimeCluster& tc, StrawHitCollection const& shcol);
      // same, taking a HelixSeed.  This uses the pitch to make a more
      // sophisticated z correction
      void updateT0(HelixSeed& hs, StrawHitCollection const& shcol);
      // access the offsets
      double strawHitTimeOffset(double hitz) const; // z position in the tracker frame!!!
      double strawHitTimeErr() const { return _shErr; }
      double caloClusterTimeOffset(int diskId) const; // depends on which 
      double caloClusterTimeErr(int diskId) const; // depends on which 
      // calculate the t0 for a straw hit. This assumes an average drift time
      double strawHitTime(StrawHit const& sh, StrawHitPosition const& shp);
      // calculate the t0 for a calo cluster.
      double caloClusterTime(CaloCluster const& cc) const;
    private:
    // helper functions
      int _debug;
//      StrawHitFlag _useflag, _dontuseflag;// flags for which hits to use
      TrkFitDirection _fdir; // fit direction.  This is used to make z-dependent time shifts
      double _shOffset;  // average time offset for straw hits
      double _shSlope; // dt/dz for straw hits
      double _shErr;
      double _caloT0Offset[2]; // time offsets for downstream particls in the 2 disks
      double _caloT0Err[2]; // time offsets errors for downstream particls in the 2 disks
  };
}
#endif
