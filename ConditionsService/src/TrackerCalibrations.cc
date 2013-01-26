//
// Parameters for tracker calibrations.
//
// $Id: TrackerCalibrations.cc,v 1.13 2013/01/26 18:16:38 brownd Exp $
// $Author: brownd $
// $Date: 2013/01/26 18:16:38 $
//

// Mu2e include files
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "TrackerGeom/inc/Straw.hh"

namespace mu2e {

  TrackerCalibrations::TrackerCalibrations( SimpleConfig const& config ){

    // Here we should eventually interface to some database
    _tdresopar0 = config.getDouble("TDResolution_0",64.2);
    _tdresopar1 = config.getDouble("TDResolution_1",60.7);
    // simplistic placeholder for drift calibration parameters
    _vdrift = config.getDouble("DriftVelocity",0.05); // mm/ns
    _rres = config.getDouble("DriftRadiusResolution",0.1); //mm
    _distvsdeltat = config.getDouble("SignalVelocity",231.); //mm/ns
    _edepToAmpl = config.getDouble("EdepToAmpl",1.0); // mV/MeV
    _amplRes = config.getDouble("AmplRes", 0.0); //   relative
    _crossTalk = config.getDouble("Crosstalk",0.0); //   relative
  }
  
  void TrackerCalibrations::DistanceToTime(StrawIndex strawIndex,double rdrift, CLHEP::Hep3Vector const& tdir,D2T& d2t) const {
    // oversimplfied model, FIXME!!!
    // Note that negative drift radii are allowed: this is necessary to allow continuous derivatives at the wire.
    // Calling classes that require a positive time should pass abs(rdrift).
    d2t._tdrift = rdrift/_vdrift;
    d2t._tdrifterr = _rres/_vdrift;
    d2t._vdrift = _vdrift;
  }

  void TrackerCalibrations::TimeToDistance(StrawIndex strawIndex, double tdrift, CLHEP::Hep3Vector const& tdir,T2D& t2d) const {
    t2d._rdrift = tdrift*_vdrift;
    t2d._rdrifterr = _rres;
    t2d._vdrift = _vdrift;
  }


  void TrackerCalibrations::EnergyToAmplitude(StrawIndex strawIndex, double edep, E2A& e2a) const {
    // oversimplfied model, FIXME!!!
    e2a._ampl = edep/_edepToAmpl;
    e2a._amplerr = e2a._ampl*_amplRes;
  }

  void TrackerCalibrations::AmplitudeToEnergy(StrawIndex strawIndex, double ampl, A2E& a2e) const {
    // oversimplfied model, FIXME!!!
    a2e._edep = ampl*_edepToAmpl;
    a2e._edeperr = ampl*_edepToAmpl*_amplRes;

  }

  double TrackerCalibrations::TimeDivisionResolution(StrawIndex , double znorm) const {
    double reso  = _tdresopar0 + _tdresopar1 * (znorm - 0.5) * (znorm - 0.5); //resolution in mm
    return reso;

  }

  double TrackerCalibrations::SignalVelocity(StrawIndex ) const {
    return _distvsdeltat; //mm/ns
  }

  double TrackerCalibrations::TimeDiffToDistance(StrawIndex strawIndex, double deltaT) const{
    return 0.5 * SignalVelocity(strawIndex) * deltaT;
  }

  void TrackerCalibrations::StrawHitInfo(Straw const& straw, StrawHit const& strawhit, SHInfo&
  shinfo) const {
   // compute the position as being on the wire the distance specified by time division.  Note this can
// be beyond the physical wire!
    double vwire = SignalVelocity(strawhit.strawIndex());
    shinfo._tddist = TimeDiffToDistance(strawhit.strawIndex(),strawhit.dt());
    shinfo._pos = straw.getMidPoint() + shinfo._tddist*straw.getDirection();
// this time represents when the particle passed by the wire.
    double shlen = straw.getHalfLength();
    shinfo._time = strawhit.time() - (shlen-shinfo._tddist)/vwire;
// Position error along the wire is given by the time division
    shinfo._tdres = TimeDivisionResolution(straw.index(),0.5*(shlen-shinfo._tddist)/shlen);
// time resolution is due to intrinsic timing resolution and time difference resolution
    shinfo._timeres = shinfo._tdres/vwire;
  }


  double TrackerCalibrations::CrossTalk(StrawIndex strawIndex0, StrawIndex strawIndexN) const {
    // FIXME oversimplfied model
    return _crossTalk;
  }

}
