//
// Parameters for tracker calibrations.
//

// Mu2e include files
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "TrackerGeom/inc/Straw.hh"

namespace mu2e {

  TrackerCalibrations::TrackerCalibrations( SimpleConfig const& config ){
// Tracker Calibration information; these correspond to Vadim's June 2012 measurement using 90 Sr
// for the resolution at the straw center, and his 55 Fe measurement for the
// position dependence
// DNB 3/8/2016.  Updated to use measured LBL prototype chamber MIP resolution
    // Here we should eventually interface to some database
//    _tdresopar0 = config.getDouble("TDResolution_0",34.0);
//    _tdresopar1 = config.getDouble("TDResolution_1",32.9);
    _tdresopar0 = config.getDouble("TDResolution_0",43.0);
    _tdresopar1 = config.getDouble("TDResolution_1",0.0);
    // simplistic placeholder for drift calibration parameters
    _vdrift = config.getDouble("DriftVelocity",0.05); // mm/ns
    _rres_min = config.getDouble("MinDriftRadiusResolution",0.09); //mm
    _rres_max = config.getDouble("MaxDriftRadiusResolution",0.15); //mm
    _rres_rad = config.getDouble("DriftRadiusResolutionRadius",1.5); //mm
    _tdbuff = config.getDouble("TimeDivisionBuffer",80); //mm
    _distvsdeltat = config.getDouble("SignalVelocity",273.); //mm/ns
    _edepToAmpl = config.getDouble("EdepToAmpl",1.0); // mV/MeV
    _amplRes = config.getDouble("AmplRes", 0.0); //   relative
    _crossTalk = config.getDouble("Crosstalk",0.0); //   relative
    _mint0doca = config.getDouble("minT0DOCA", -0.2);
  }

  void TrackerCalibrations::DistanceToTime(StrawIndex strawIndex,double rdrift, CLHEP::Hep3Vector const& tdir,D2T& d2t) const {
    // oversimplfied model, FIXME!!!
    // Note that negative drift radii are allowed: this is necessary to allow continuous derivatives at the wire.
    // Calling classes that require a positive time should pass abs(rdrift).
    d2t._tdrift = rdrift/_vdrift;
    double rres = _rres_min;
    if(rdrift<_rres_rad){
      rres = _rres_min+_rres_max*(_rres_rad-rdrift)/_rres_rad;
    }
    d2t._tdrifterr = rres/_vdrift;
    d2t._vdrift = _vdrift;
  }

  void TrackerCalibrations::TimeToDistance(StrawIndex strawIndex, double tdrift, CLHEP::Hep3Vector const& tdir,T2D& t2d) const {
    t2d._rdrift = tdrift*_vdrift;
    double rres = _rres_min;
    if(t2d._rdrift<_rres_rad){
      rres = _rres_min+_rres_max*(_rres_rad-t2d._rdrift)/_rres_rad;
    }
    t2d._rdrifterr = rres;
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
  // this resolution function assumes all straws lengths have the same relative resolution FIXME!!
    double reso  = _tdresopar0 + _tdresopar1 * (znorm - 0.5) * (znorm - 0.5); //resolution in mm
    return reso;

  }

  double TrackerCalibrations::SignalVelocity(StrawIndex ) const {
    return _distvsdeltat; //mm/ns
  }

  double TrackerCalibrations::TimeDiffToDistance(StrawIndex strawIndex, double deltaT) const{
    return 0.5 * SignalVelocity(strawIndex) * deltaT;
  }

  double TrackerCalibrations::CrossTalk(StrawIndex strawIndex0, StrawIndex strawIndexN) const {
    // FIXME oversimplfied model
    return _crossTalk;
  }


  double TrackerCalibrations::Mint0doca() const{
    return _mint0doca;
  }
}
