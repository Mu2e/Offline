//
//  Persistent representation of a TrkStrawHit, used in the
//  persistent representation of the BTrk Kalman Fit
//  Original author: Dave Brown (LBNL) 31 Aug 2016
//
#ifndef RecoDataProducts_TrkStrawHitSeed_HH
#define RecoDataProducts_TrkStrawHitSeed_HH
#include "BTrk/TrkBase/TrkT0.hh"
#include "RecoDataProducts/inc/HitT0.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "DataProducts/inc/StrawId.hh"
#include <Rtypes.h>
namespace mu2e {
  struct TrkStrawHitSeed {
    TrkStrawHitSeed() : _index(0), _trklen(0), _hitlen(0), _rdrift(0), _wdoca(0), _rerr(0), _ambig(0) {}
    // construct from the information
  TrkStrawHitSeed(StrawHitIndex index, StrawId const& strawid, TrkT0 const& t0, Float_t trklen, Float_t hitlen, Float_t rdrift,
      Float_t wdoca, Int_t ambig, Float_t rerr, StrawHitFlag const& flag) :
    _index(index), _sid(strawid), _t0(t0), _trklen(trklen),
    _hitlen(hitlen), _rdrift(rdrift), _wdoca(wdoca), _rerr(rerr), _ambig(ambig), _flag(flag)  {}
// accessors
    StrawHitIndex	index() const { return _index; }
    StrawId const&	strawId() const { return _sid; }
    TrkT0 const&  t0() const { return _t0; }
    Float_t	trkLen() const { return _trklen; }
    Float_t	hitLen() const { return _hitlen; }
    Float_t	driftRad() const { return _rdrift; }
    Float_t	wireDOCA() const { return _wdoca; }
    Float_t	radialErr() const { return _rerr; }
    Int_t	ambig() const { return _ambig; }
    StrawHitFlag const& flag() const { return _flag; }
    //
    StrawHitIndex	    _index;	  // index into the primary StrawHit (and StrawHitPosition, StrawHitFlag, ..) collection
    StrawId	    _sid;	  // which straw has the hit
    TrkT0	    _t0;	  // time origin for this hit
    Float_t	    _trklen;	  // Length from the nominal track start to the POCA of this hit
    Float_t	    _hitlen;	  // Length from the straw center to the POCA of this hit
    Float_t	    _rdrift;	  // drift radius for this hit
    Float_t	    _wdoca;	  // DOCA from the track to the wire, signed by the angular momentum WRT the wire
    Float_t	    _rerr;	  // intrinsic radial error
    Int_t	    _ambig;	  // LR ambiguity assigned to this hit
    StrawHitFlag    _flag;	  // flag describing the status of this hit (active, ....)
  };
}
#endif
