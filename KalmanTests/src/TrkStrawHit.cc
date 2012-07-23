//
// BaBar hit object corresponding to a single straw hit
//
// $Id: TrkStrawHit.cc,v 1.20 2012/07/23 17:52:27 brownd Exp $
// $Author: brownd $ 
// $Date: 2012/07/23 17:52:27 $
//
// Original author David Brown, LBNL
//
#include "BaBar/BaBar.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/TrkPoca.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "TrkBase/TrkDetElemId.hh"
#include "TrkBase/TrkRep.hh"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include <algorithm>

using namespace std;

namespace mu2e
{
// convert speed of light in mm/nsec
  const double TrkStrawHit::_vlight = CLHEP::c_light;
/// Material information, BaBar style
  MatDBInfo* TrkStrawHit::matdbinfo(){
    static MatDBInfo mat;
    return &mat;
  }
  DetStrawHitType* TrkStrawHit::wtype(){
    static DetStrawHitType instance(matdbinfo(),"straw-wall");
    return &instance;
  }
  DetStrawHitType* TrkStrawHit::gtype(){
    static DetStrawHitType instance(matdbinfo(),"straw-gas");
    return &instance;
  }

  TrkStrawHit::TrkStrawHit(const StrawHit& strawhit, const Straw& straw, unsigned istraw,
    const TrkT0& trkt0,double flt0,double fltlen,double exterr,double maxdriftpull) :
    TrkHitOnTrk(1e-5),
    _strawhit(strawhit),
    _straw(straw),
    _istraw(istraw),
    _exterr(exterr),
    _penerr(0.0),
    _toterr(0.0),
    _iamb(0),
    _ambigupdate(false),
    _maxdriftpull(maxdriftpull),
    _welem(wtype(),this),
    _gelem(gtype(),this)
  {
// is there an efficiency issue fetching the calibration object for every hit???
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    tcal->StrawHitInfo(strawhit,_wpos,_wtime,_tddist_err,_wtime_err);
    CLHEP::Hep3Vector const& wiredir = _straw.getDirection();
// get time division and drift information for this straw hit relative to the wire center
    _tddist = tcal->TimeDiffToDistance(_straw.index(),_strawhit.dt());
    CLHEP::Hep3Vector const& mid = _straw.getMidPoint();
// the hit trajectory is defined as a line segment directed along the wire direction starting from the wire center
    _hittraj = new TrkLineTraj(HepPoint(mid.x(),mid.y(),mid.z()),wiredir,_tddist-_tddist_err,_tddist+_tddist_err);
    setHitLen(_tddist);
    setFltLen(fltlen);
// compute initial hit t0 and drift
    updateT0(trkt0,flt0);
//    std::cout << "creating TrkStrawHit " << this << std::endl;
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
      _exterr(other._exterr),
      _penerr(other._penerr),
      _toterr(other._toterr),
      _iamb(other._iamb),
      _ambigupdate(false),
      _t2d(other._t2d),
      _tddist(other._tddist),
      _tddist_err(other._tddist_err),
      _wtime(other._wtime),
      _wtime_err(other._wtime_err),
      _maxdriftpull(other._maxdriftpull),
      _welem(wtype(),this),
      _gelem(gtype(),this)
  {
//    std::cout << "creating TrkStrawHit copy " << this << std::endl;
  }
  
  TrkStrawHit::~TrkStrawHit(){
// delete the hit
//    std::cout << "deleting hit " << _theHit << std::endl;
    delete _hittraj;
// ugly trick to keep the base class from trying to delete _TrkDummyHit
    _parentRep=0;  
//    std::cout << "deleted TrkStrawHit " << this << std::endl;
  }

  TrkStrawHit*
  TrkStrawHit::clone(TrkRep* parentRep, const TrkDifTraj* trkTraj) const {
    return new TrkStrawHit(*this, parentRep);
  }
  
  void
  TrkStrawHit::updateT0(const TrkT0& trkt0,double t0flt){
// correct for track propagation and signal propagation
    double dflt = fltLen() - t0flt;
// the following assumes beta=1, this should come from the fit, FIXME!!!
    double tprop = dflt/_vlight;  
    // hitt0 is the time the particle passed closest to the wire.  It does not include
// drift time or electronics signal propagation time.
    _hitt0 = trkt0._t0 + tprop;
// add a term that scales linearly with the flight distance from t0flt.  This should be taken
// out by the Kalman fit, but it needs additional derivative terms.  FIXME!!!!
    _hitt0_err = trkt0._t0err;
  }
  
  double
  TrkStrawHit::time() const {
    return strawHit().time();
  }
  
  void
  TrkStrawHit::updateDrift() {
// compute the electronics propagation time.  The convention is that the hit time is measured at the
// FAR END of the wire, as signed by the wire direction.
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    double vwire = tcal->SignalVelocity(straw().index());
    double teprop;
    if(_poca != 0 && _poca->status().success()){
      teprop = (straw().getHalfLength()-hitLen())/vwire;
    } else {
// if we're missing poca information, use time division instead
      teprop = (straw().getHalfLength()-_tddist)/vwire;
    }
    double tdrift = strawHit().time() - _hitt0 - teprop;
// find the track direction at this hit
    CLHEP::Hep3Vector tdir = getParentRep()->traj().direction(fltLen());
// convert time to distance.  This computes the intrinsic drift radius error as well
    tcal->TimeToDistance(straw().index(),tdrift,tdir,_t2d);
// Propogate error in t0, using local drift velocity
    double rt0err = _hitt0_err*_t2d._vdrift;
// deal with ambiguity updating
    if(_ambigupdate) {
      int iamb = _poca->doca() > 0 ? 1 : -1;
      setAmbig(iamb);
// pply a penalty if the drift distance is too small
//      if(fabs(_t2d._rdrift) < _
    }
    // total hit error is the sum of all
    _toterr = sqrt(_t2d._rdrifterr*_t2d._rdrifterr + rt0err*rt0err + _exterr*_exterr + _penerr*_penerr);
// If the hit is wildly away from the track , disable it
    double rstraw = _straw.getRadius(); 
    if( _t2d._rdrift - rstraw > _maxdriftpull*_toterr ||
      _t2d._rdrift < -_maxdriftpull*_toterr){
        setUsability(-10);
        setActivity(false);
    } else {
// otherwise restrict to a physical range
      if (_t2d._rdrift < 0.0){
	_t2d._rdrift = 0.0;
      } else if( _t2d._rdrift > rstraw){
	_t2d._rdrift = rstraw;
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
  TrkStrawHit::updateMeasurement(const TrkDifTraj* traj) {
    TrkErrCode status(TrkErrCode::fail);
// find POCA to the wire
    updatePoca(traj);
   if(_poca != 0 && _poca->status().success()) {
      status = _poca->status();
// update the drift distance using this traj direction
      updateDrift();
 // sign drift distance by ambiguity.  Note that an ambiguity of 0 means to ignore the drift
      double residual = _poca->doca() - _t2d._rdrift*_iamb;
      setHitResid(residual);
      setHitRms(_toterr);
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
      hpos = _wpos + pdir*_t2d._rdrift*_iamb;
    } else {
      hpos = _wpos;
    }
  }

// compute the pathlength through one wall of the straw, given the drift distance and straw geometry
  double
  TrkStrawHit::wallPath(double pdist,Hep3Vector const& tdir) const {
    double thick = straw().getThickness();
    double radius = straw().getRadius();
    if(pdist>=radius)
      pdist = 0.96*radius;
    double wallpath =  (sqrt( (radius+thick+pdist)*(radius+thick-pdist) ) -
      sqrt( (radius+pdist)*(radius-pdist) ));
  // scale for the other dimension
    double cost = tdir.dot(_straw.getDirection());
    if(fabs(cost)<0.999)
      wallpath /= sqrt( (1.0-cost)*(1.0+cost) );
    else
      wallpath = radius;
// use half-length as maximum length
    wallpath = std::min(wallpath,radius);
// test for NAN	
//    if(wallpath != wallpath){
///      std::cout << "non wall" << std::endl;
//   }
    return wallpath;
  }
  
  // compute the pathlength through half the gas , given the drift distance and straw geometry
  double
  TrkStrawHit::gasPath(double pdist,Hep3Vector const& tdir) const {
    double radius = straw().getRadius();
    double hlen = straw().getHalfLength();
    if(pdist>=radius)
      pdist = 0.96*radius;
    double gaspath = sqrt( (radius+pdist)*(radius-pdist) );
// scale for the other dimension
    double cost = tdir.dot(_straw.getDirection());
    if(fabs(cost)<0.999)
      gaspath /= sqrt( (1.0-cost)*(1.0+cost) );
    else
      gaspath = hlen;
// use half-length as maximum length
    gaspath = std::min(gaspath,hlen);
//NAN test
//    if(gaspath != gaspath){
//      std::cout << "nan gas" << endl;
//    }
    return gaspath;
  }

} 
