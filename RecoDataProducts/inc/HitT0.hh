//
// Class representing the t0 of a given hit: this is the time the physical particle came
// closest to the wire, and is the offset for measuring the drift time.
//  Original Author: Dave Brown (LBNL) 21 Aug. 2016
//
#ifndef RecoDataProducts_HitT0_HH
#define RecoDataProducts_HitT0_HH
#include <Rtypes.h>
#include "Offline/BTrkLegacy/inc/TrkT0.hh"
namespace mu2e {
  struct HitT0 {
    Float_t _t0; // t0 value
    Float_t _t0err; // error on t0
    Float_t t0() const { return _t0; }
    Float_t t0Err() const { return _t0err; }
    HitT0() : _t0(0.0), _t0err(0.0) {}
    // conversion with TrkT0
    HitT0(TrkT0 const& trkt0) : _t0(trkt0.t0()), _t0err(trkt0.t0Err()) {}
    HitT0(double t0, double t0err) : _t0(t0), _t0err(t0err) {}
    HitT0& operator =(TrkT0 const& trkt0) { _t0 = trkt0.t0(); _t0err = trkt0.t0Err(); return *this; }
    operator TrkT0() const { return TrkT0(_t0,_t0err); }
  };
}
#endif
