//
// BaBar hit object corresponding to a single straw hit
//
// $Id: TrkStrawHit.cc,v 1.8 2011/07/13 20:44:27 mu2ecvs Exp $
// $Author: mu2ecvs $ 
// $Date: 2011/07/13 20:44:27 $
//
// Original author David Brown, LBNL
//
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/TrkPoca.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "TrkBase/TrkDetElemId.hh"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

using namespace std;

namespace mu2e
{
  
  // drift velocity is a parameter of the makeStrawHits module; I have no access to that, so I hardcode
  // this here.  I'm also hard-coding the wire signal propagation velocity.  FIXME!!!!!
  double TrkStrawHit::_vdrift = 0.05; // mm per nanosecond
  double TrkStrawHit::_maxdriftpull = 4.0; // disable hits with unphysical drift distance pulls beyond this cut

  TrkDummyHit::TrkDummyHit(TrkEnums::TrkViewInfo v, int id, TrkDetElemId::systemIndex sys)
      : _view(v), _eid(id,sys)
  {}
  TrkDummyHit::TrkDummyHit(const TrkDummyHit& other)
      : _view(other._view), _eid(other._eid)
  {}
  TrkDummyHit::~TrkDummyHit() {}
  TrkDummyHit* TrkDummyHit::clone() const { return new TrkDummyHit(*this);}
  
  
  TrkStrawHit::TrkStrawHit(const StrawHit& strawhit, const Straw& straw, unsigned istraw, double hitt0, double hitt0_err,double herr) :
    TrkHitOnTrk(new TrkDummyHit(TrkEnums::xyView,strawhit.strawIndex().asInt(),TrkDetElemId::null),1e-5),
    _strawhit(strawhit),
    _straw(straw),
    _istraw(istraw),
    _hitt0(hitt0),
    _hitt0_err(hitt0_err),
    _herr(herr),
    _iamb(0)
  {
// is there an efficiency issue fetching the calibration object for every hit???
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    tcal->StrawHitInfo(strawhit,_wpos,_wtime,_tddist_err,_wtime_err);
    CLHEP::Hep3Vector const& wiredir = _straw.getDirection();
  // get time division and drift information for this straw hit relative to the wire center
    _tddist = tcal->TimeDiffToDistance(_straw.index(),_strawhit.dt());
    CLHEP::Hep3Vector const& mid = _straw.getMidPoint();
    _hittraj = new TrkLineTraj(HepPoint(mid.x(),mid.y(),mid.z()),wiredir,_tddist-_tddist_err,_tddist+_tddist_err);
    setHitLen(_tddist);
  // compute drift parameters, if the initial t0 error was positive
    if(_hitt0_err>0.0)updateDrift();
  }
  
  TrkStrawHit::TrkStrawHit(const TrkStrawHit& other, TrkRep* rep) :
      TrkHitOnTrk(other, rep, 0),
      _strawhit(other._strawhit),
      _straw(other._straw),
      _istraw(other._istraw),
      _hittraj(other._hittraj->clone()),
      _wpos(other._wpos),
      _hitt0(other._hitt0),
      _hitt0_err(other._hitt0_err),
      _herr(other._herr),
      _iamb(other._iamb),
      _rdrift(other._rdrift),
      _rdrift_err(other._rdrift_err),
      _tddist(other._tddist),
      _tddist_err(other._tddist_err)
  {
  }
  
  TrkStrawHit::~TrkStrawHit(){
//    delete _hit;
    delete _hittraj;
// ugly trick to keep the base class from trying to delete _TrkDummyHit
    _parentRep=0;  
  }

  TrkStrawHit*
  TrkStrawHit::clone(TrkRep* parentRep, const TrkDifTraj* trkTraj) const {
    return new TrkStrawHit(*this, parentRep);
  }
  
  void
  TrkStrawHit::updateT0(double hitt0,double hitt0_err){
    _hitt0 = hitt0;
    _hitt0_err = hitt0_err;
    updateDrift();
  }
  
  double
  TrkStrawHit::time() const {
    return _wtime;
  }
  
  void
  TrkStrawHit::updateDrift() {
// assume wire propagation velocity is speed of light: FIXME!!!
    double tdrift = time() - _hitt0;
    _rdrift = tdrift*_vdrift;
  // radial error is combination of time error and intrinsic error.  This doesn't account for correlation between hits
    double time_err = sqrt(_hitt0_err*_hitt0_err + _wtime_err*_wtime_err);
    double rt0err = time_err*_vdrift;
    _rdrift_err = sqrt(rt0err*rt0err + _herr*_herr);
    if(_rdrift > _straw.getRadius()){
  // unphysical condition: decide if this hit should get disabled, or just brought
  // back into the physical range.  Parameters here should be adjustable FIXME!!!
      double dr = _rdrift - _straw.getRadius();
      if( dr < _maxdriftpull*_rdrift_err){
//        _rdrift = _straw.getRadius();
      } else {
        setUsability(-10);
        setActivity(false);
      }
    } else if (_rdrift < 0.0){
      if( fabs(_rdrift) < _maxdriftpull*_rdrift_err){
//        _rdrift = 0.0;
      } else {
        setUsability(-10);
        setActivity(false);
      }
    }
  }
  
  void 
  TrkStrawHit::setAmbig(int newambig){
//    if(newambig != _iamb && _iamb != 0)
//      std::cout << "changing hit ambiguity " << std::endl;
    if(newambig > 0)
      _iamb = 1;
    else if(newambig < 0)
      _iamb = -1;
    else
      _iamb = 0;
  }

  TrkErrCode
  TrkStrawHit::updateMeasurement(const TrkDifTraj* traj, bool maintainAmbiguity) {
    TrkErrCode status(TrkErrCode::fail);
// find POCA to the wire
    updatePoca(traj, maintainAmbiguity);
    if(_poca != 0 && _poca->status().success()) {
      status = _poca->status();
// set the ambiguity if it was never set, or if allowed, check both ambiguities and set as necessary
      if(_iamb == 0 || !maintainAmbiguity){
// reset ambiguity only if the difference is significant.  This avoids frothing during the fit iterations
        if(_rdrift > _herr ||  _iamb == 0){
          int newamb = _poca->doca() > 0 ? 1 : -1;
          setAmbig(newamb);
        }
      }
// sign drift distance by ambiguity
      double residual = _poca->doca() - _rdrift*_iamb;
      setHitResid(residual);
      setHitRms(_rdrift_err);
    } else {
      cout << "TrkStrawHit:: updateMeasurement() failed" << endl;
      setHitResid(999999);
      setHitRms(999999);
      setUsability(0);
    }
    return status;
  }
  
  void
  TrkStrawHit::hitPosition(CLHEP::Hep3Vector& hpos) const{
    if(_poca != 0 && _poca->status().success()){
      CLHEP::Hep3Vector pdir = (trkTraj()->position(fltLen()) - hitTraj()->position(hitLen())).unit();
      hpos = _wpos + pdir*_rdrift*_iamb;
    } else {
      hpos = _wpos;
    }
  }

// compute the pathlength through one wall of the straw, given the drift distance and straw geometry
  double
  TrkStrawHit::wallPath() const {
    double thick = straw().getThickness();
    double radius = straw().getRadius();
    double drift = min(radius,driftRadius());
    double wallpath =  sqrt( (radius+thick+drift)*(radius+thick-drift) ) -
      sqrt( (radius+drift)*(radius-drift) );
    return wallpath;
  }
  
  // compute the pathlength through one wall of the straw, given the drift distance and straw geometry
  double
  TrkStrawHit::gasPath() const {
    double radius = straw().getRadius();
    double drift = min(radius,driftRadius());
    double gaspath = sqrt( (radius+drift)*(radius-drift) );
    return gaspath;
  }

} 
