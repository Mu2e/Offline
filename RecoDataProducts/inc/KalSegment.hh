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
#include "RecoDataProducts/inc/HitT0.hh"
#include "KinKal/General/ParticleState.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
#include "KinKal/Trajectory/CentralHelix.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "CLHEP/Matrix/Vector.h"

namespace mu2e {
  struct KalSegment {
    KalSegment() : _tmin(0.0), _tmax(-1.0), _dflt(0.0) {}
    KalSegment(KinKal::TimeRange const& trange) : _tmin(trange.begin()), _tmax(trange.end()){}
    // construct from a kinematic trajectory
    template <class KTRAJ> explicit KalSegment(KTRAJ const& ktraj, double time, double dflt=0.0) :
      _tmin(ktraj.range().begin()), _tmax(ktraj.range().end()),
      _bnom(ktraj.bnom()),
      _pstate(ktraj.stateEstimate(time)) , _dflt(dflt) {}
    // accessors
    KinKal::ParticleStateEstimate const& state() const { return _pstate; }
    // provide a few call-down functions
    double mom() const { return _pstate.momentum(); }
    double momerr() const { return sqrt(_pstate.momentumVariance()); }
    KinKal::VEC3 momentum3() const { return _pstate.momentum3(); }
    KinKal::VEC3 position3() const { return _pstate.position3(); }
 //   void mom(float fltlen, XYZVec& momvec) const { helix().direction(fltlen,momvec); momvec *= mom(); } // momentum as a function of local flightlength
    // convert content to a LoopHelix
    KinKal::LoopHelix loopHelix() const { return KinKal::LoopHelix(_pstate, KKbnom(),KinKal::TimeRange(tmin(),tmax())); }
    // convert to a CentralHelix
    KinKal::CentralHelix centralHelix() const { return KinKal::CentralHelix(_pstate, KKbnom(),KinKal::TimeRange(tmin(),tmax())); }
    KinKal::KinematicLine kinematicLine() const { return KinKal::KinematicLine(_pstate, KKbnom(),KinKal::TimeRange(tmin(),tmax())); }
    Float_t tmin() const { return _tmin; }
    Float_t tmax() const { return _tmax; }
    auto tref() const { return _pstate.time(); }
    // t0 = time (and error) for when particle goes through z=0; 
    HitT0 t0() const;
    XYZVec const& bnom() const { return _bnom; }
    KinKal::VEC3 KKbnom() const { return KinKal::VEC3(_bnom); }
    Float_t _tmin, _tmax; // time range
// main payload is the particle state estimate.  this includes all the kinematic information to
// interpret as anything else.  BField is needed to interpret geometrically
    XYZVec _bnom; // Bfield associated with this segment, needed to reconstitute helix
    KinKal::ParticleStateEstimate _pstate; // particle state at this sample
 // the following are deprecated legacy functions specific to the BTrk fit and will go away eventually
    HelixVal helix() const;
    HelixCov covar() const;
    void mom(double flt, XYZVec& momvec) const;
    double fmin() const { return timeToFlt(_tmin); } // local 3D flight range
    double fmax() const { return timeToFlt(_tmax); }
    Float_t _dflt;    
    double localFlt(float globalflt) const { return globalflt + _dflt; }
    double globalFlt(float localflt) const { return localflt - _dflt; }
    double fltToTime(double flt) const; // local flight
    double timeToFlt(double time) const; // local flight
  };
}
#endif

