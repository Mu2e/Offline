//
//  Persistent representation of the BTrk Kalman filter fit (KalRep)
//  Oiriginal author: Dave Brown (LBNL) 31 Aug 2016
//
#ifndef RecoDataProducts_KalSeed_HH
#define RecoDataProducts_KalSeed_HH
// mu2e
#include "RecoDataProducts/inc/TrkFitDirection.hh"
#include "RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "RecoDataProducts/inc/KalSegment.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
// BTrk
#include "BTrk/TrkBase/TrkParticle.hh"
// Root
#include <Rtypes.h>
// C++
#include <vector>
namespace mu2e {
  struct KalSeed {
    KalSeed() :  _flt0(0) {}
    KalSeed(TrkParticle tpart,TrkFitDirection fdir,TrkT0 const& t0, double flt0, TrkFitFlag const& status) :
      _tpart(tpart), _fdir(fdir), _t0(t0), _flt0(static_cast<Float_t>(flt0)), _status(status) {}


    TrkParticle const& particle() const { return _tpart; }
    TrkFitDirection const& fitDirection() const { return _fdir; }
    std::vector<TrkStrawHitSeed> const& hits() const { return _hits;}
    std::vector<KalSegment> const& segments() const { return _segments; }
    TrkFitFlag const& status() const { return _status; }
    Float_t flt0() const { return _flt0; }
    TrkT0 const& t0() const { return _t0; }

    TrkParticle			    _tpart; // particle assumed for this fit
    TrkFitDirection	      	    _fdir; // direction in which this particle was fit
    TrkT0			    _t0; // track t0
    Float_t			    _flt0; // flight distance where the track crosses the tracker midplane (z=0)
    std::vector<TrkStrawHitSeed>    _hits; // hit seeds for all the hits used in this fit
    std::vector<KalSegment>	    _segments; // segments of the Kalman filter fit result
    TrkFitFlag			    _status; // status of this fit
    // eventually add fit quality information FIXME!
  };
}
#endif
