//
//  Persistent representation of a TrkStrawHit, used in the
//  persistent representation of the BTrk Kalman Fit
//  Original author: Dave Brown (LBNL) 31 Aug 2016
//
#ifndef RecoDataProducts_TrkStrawHitSeed_HH
#define RecoDataProducts_TrkStrawHitSeed_HH
#include "RecoDataProducts/inc/HitT0.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "DataProducts/inc/StrawId.hh"
#include <Rtypes.h>
#include <functional>
namespace mu2e {
  struct TrkStrawHitSeed {
    TrkStrawHitSeed() : _index(0), _trklen(0), _hitlen(0), _rdrift(0), _dtime(0), _stime(0), _htime(0),
      _wdoca(0), _rerr(0), _ambig(0), _edep(0), _wdist(0), _werr(0) {}
    // construct from the information
    TrkStrawHitSeed(StrawHitIndex index, HitT0 const& t0, Float_t trklen, Float_t hitlen, Float_t rdrift,
	Float_t stime,
	Float_t wdoca, Int_t ambig, Float_t rerr, StrawHitFlag const& flag, ComboHit const& chit) :
      _index(index), _sid(chit.strawId()), _t0(t0), _trklen(trklen),
      _hitlen(hitlen), _rdrift(rdrift),
      _dtime(chit.driftTime()), _stime(stime),_htime(chit.time()),
      _wdoca(wdoca), _rerr(rerr), _ambig(ambig), 
      _edep(chit.energyDep()),_wdist(chit.wireDist()), _werr(chit.wireRes()), _end(chit.driftEnd()), 
      _flag(flag)  {}
    TrkStrawHitSeed(StrawHitIndex index, TrkT0 const& t0, Float_t trklen, Float_t hitlen, Float_t rdrift,
	Float_t stime,
	Float_t wdoca, Int_t ambig, Float_t rerr, StrawHitFlag const& flag, ComboHit const& chit) :
      _index(index), _sid(chit.strawId()), _t0(t0), _trklen(trklen),
      _hitlen(hitlen), _rdrift(rdrift),
      _dtime(chit.driftTime()), _stime(stime),_htime(chit.time()),
      _wdoca(wdoca), _rerr(rerr), _ambig(ambig), 
      _edep(chit.energyDep()),_wdist(chit.wireDist()), _werr(chit.wireRes()), _end(chit.driftEnd()), 
      _flag(flag)  {}
    // accessors
    StrawHitIndex index() const { return _index; }
    StrawId const&	strawId() const { return _sid; }
    HitT0 const&  t0() const { return _t0; }
    Float_t	trkLen() const { return _trklen; }
    Float_t	hitLen() const { return _hitlen; }
    Float_t	driftRadius() const { return _rdrift; }
    Float_t	signalTime() const { return _stime; }
    Float_t	hitTime() const { return _htime; }
    Float_t	TOTDriftTime() const { return _dtime; }
    Float_t	wireDOCA() const { return _wdoca; }
    Float_t	radialErr() const { return _rerr; }
    Float_t     wireDist() const { return _wdist; }
    Float_t     wireRes() const { return _werr; }
    Float_t     energyDep() const { return _edep; }
    StrawEnd    driftEnd() const { return _end; }
    Int_t	ambig() const { return _ambig; }
    StrawHitFlag const& flag() const { return _flag; }
    //
    StrawHitIndex   _index;       // index to the original straw (Combo) hit, and (for MC) MCDigi
    StrawId	    _sid;	  // which straw has the hit
    HitT0	    _t0;	  // time origin for this hit = track t0 + particle propagation time to this straw
    Float_t	    _trklen;	  // Length from the nominal track start to the POCA of this hit
    Float_t	    _hitlen;	  // Length from the straw center to the POCA of this hit
    Float_t	    _rdrift;	  // drift radius for this hit
    Float_t	    _dtime;	  // drift time from TOT for this hit
    Float_t	    _stime;	  // signal propagation time for this hit, to the nearest end
    Float_t	    _htime;	  // measured time for this hit 
    Float_t	    _wdoca;	  // DOCA from the track to the wire, signed by the angular momentum WRT the wire
    Float_t	    _rerr;	  // intrinsic radial error
    Int_t	    _ambig;	  // LR ambiguity assigned to this hit
    Float_t         _edep;        // reco energy deposition
    Float_t         _wdist;       // wire distance
    Float_t         _werr;        // wire distance error estimate
    StrawEnd	    _end;         // straw end used for hit time measurement
    StrawHitFlag    _flag;	  // flag describing the status of this hit (active, ....)
  };
  // binary functor to sort TrkStrawHits by StrawHit index
  struct indexcompseed : public std::binary_function<TrkStrawHitSeed,TrkStrawHitSeed, bool> {
    bool operator()(const TrkStrawHitSeed& x,const TrkStrawHitSeed& y) { return x.index() < y.index(); }
  };
}
#endif
