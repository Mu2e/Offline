//
//  Class representing a segment of the Kalman filter track fit.  The segment contains part of the fit
//  result, accurate and unchanging over the specified range of the particle trajectory.
//  Segments can be strung together to represent the full Kalman fit, or they can be sampled
//  to provide valid fit results in subset of positions.
//  Original Author: Dave Brown (LBNL) 31 Aug. 2016
//
#ifndef RecoDataProducts_KalSegment_HH
#define RecoDataProducts_KalSegment_HH
#include <Rtypes.h>
#include "RecoDataProducts/inc/HelixVal.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "KinKal/General/ParticleState.hh"
#include "KinKal/Trajectory/CentralHelix.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "CLHEP/Matrix/Vector.h"

namespace mu2e {
  struct KalSegment {
    KalSegment() : _tmin(0.0), _tmax(-1.0), _fmin(0.), _fmax(-1.0), _dflt(0.0) {}
    // provide a few call-down functions
    double mom() const { return _pstate.momentum(); }
    double momerr() const { return sqrt(_pstate.momentumVariance()); }
    KinKal::VEC3 momentum3() const { return _pstate.momentum3(); }
    KinKal::VEC3 position3() const { return _pstate.position3(); }
 //   void mom(float fltlen, XYZVec& momvec) const { helix().direction(fltlen,momvec); momvec *= mom(); } // momentum as a function of local flightlength
    // conver the content as a LoopHelix
    KinKal::LoopHelix loopHelix() const { return KinKal::LoopHelix(_pstate, bnom(),KinKal::TimeRange(tmin(),tmax())); }
    // same, as a CentralHelix
    KinKal::CentralHelix centralHelix() const { return KinKal::CentralHelix(_pstate, bnom(),KinKal::TimeRange(tmin(),tmax())); }
    Float_t tmin() const { return _tmin; }
    Float_t tmax() const { return _tmax; }
    KinKal::VEC3 bnom() const { return KinKal::VEC3(_bnom.X(), _bnom.Y(), _bnom.Z()); }
    XYZVec _bnom; // Bfield associated with this segment, needed to reconstitute helix
    Float_t _tmin, _tmax; // time range
// main payload is the particle state estimate.  this includes all the kinematic information to
// interpret as anything else.  BField is needed to interpret geometrically
    KinKal::ParticleStateEstimate _pstate; // particle state at this sample
 // the following are deprecated legacy functions specific to the BTrk fit and will go away eventually
    HelixVal helix() const;
    HelixCov covar() const;
    void mom(double flt, XYZVec& momvec) const;
    Float_t fmin() const { return _fmin; }
    Float_t fmax() const { return _fmax; }
    float localFlt(float globalflt) const { return globalflt + _dflt; }
    float globalFlt(float localflt) const { return localflt - _dflt; }
    double fltToTime(float flt) const { return _tmin + (_tmax-_tmin)*(flt-_fmin)/(_fmax - _fmin); }
    // legacy payload for BTrk
    Float_t _fmin, _fmax; // legacy for BTrk
    Float_t _dflt; // legacy
  };
}
#endif

