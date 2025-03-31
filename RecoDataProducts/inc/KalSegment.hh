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
#include "Offline/RecoDataProducts/inc/HelixVal.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/RecoDataProducts/inc/HitT0.hh"
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"
#include "KinKal/General/ParticleStateEstimate.hh"
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
    KinKal::VEC3 velocity() const { return _pstate.velocity(); }
    KinKal::VEC3 position3() const { return _pstate.position3(); }
    // convert content to a LoopHelix
    KinKal::LoopHelix loopHelix() const { return KinKal::LoopHelix(_pstate, bnom(),timeRange()); }
    // convert to a CentralHelix
    KinKal::CentralHelix centralHelix() const { return KinKal::CentralHelix(_pstate, bnom(),timeRange()); }
    // convert to a KinematicLine
    KinKal::KinematicLine kinematicLine() const { return KinKal::KinematicLine(_pstate, bnom(),timeRange()); }
    double const& tmin() const { return _tmin; }
    double const& tmax() const { return _tmax; }
    KinKal::TimeRange timeRange() const;
    auto tref() const { return _pstate.time(); }// note this is NOT t0; it is the reference time at which the trajectory was sampled
    double t0Val(TrkFitFlag const& flag = TrkFitFlag(TrkFitFlag::KKLoopHelix)) const;// this will give the local segment's t0, interpreted for the given trajectory type,  which is not necessarily the same as the track t0
    KinKal::VEC3 const& bnom() const { return _bnom; }
    double _tmin, _tmax; // time range
// main payload is the particle state estimate.  this includes all the kinematic information to
// interpret as anything else.  BField is needed to interpret geometrically
    KinKal::VEC3 _bnom; // Bfield associated with this segment, needed to reconstitute helix. Needs to be double precision to avoid
    // roundoff when reconstituting parameers
    KinKal::ParticleStateEstimate _pstate; // particle state at this sample
 // the following are deprecated legacy functions specific to the BTrk, these should be removed FIXME
    HitT0 t0() const;
    HelixVal helix() const;
    HelixCov covar() const;
    void mom(double flt, XYZVectorF& momvec) const;
    double fmin() const { return timeToFlt(_tmin); } // local 3D flight range
    double fmax() const { return timeToFlt(_tmax); }
    double _dflt;
    double localFlt(double globalflt) const { return globalflt + _dflt; }
    double globalFlt(double localflt) const { return localflt - _dflt; }
    double fltToTime(double flt) const; // local flight
    double timeToFlt(double time) const; // local flight
  };
}
#endif

