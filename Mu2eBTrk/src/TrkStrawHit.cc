//
// BaBar hit object corresponding to a single straw hit
//
// $Id: TrkStrawHit.cc,v 1.25 2014/08/22 16:10:41 tassiell Exp $
// $Author: tassiell $ 
// $Date: 2014/08/22 16:10:41 $
//
// Original author David Brown, LBNL
//
#include "Mu2eBTrk/inc/TrkStrawHit.hh"
// BTrk
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/TrkDifTraj.hh"
#include "BTrk/TrkBase/TrkDetElemId.hh"
#include "BTrk/TrkBase/TrkRep.hh"
//
#include "CLHEP/Vector/ThreeVector.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include <algorithm>

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e
{
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
    const TrkT0& trkt0,double fltlen,double exterr,double maxdriftpull) :
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
    SHInfo shinfo;
    tcal->StrawHitInfo(straw,strawhit,shinfo);
    _wpos = shinfo._pos;
    _tddist = shinfo._tddist;
    _tddist_err = shinfo._tdres;
    Hep3Vector const& wiredir = _straw.getDirection();
    Hep3Vector const& mid = _straw.getMidPoint();
// the hit trajectory is defined as a line segment directed along the wire direction starting from the wire center
    _hittraj = new TrkLineTraj(HepPoint(mid.x(),mid.y(),mid.z()),wiredir,_tddist-_tddist_err,_tddist+_tddist_err);
    setHitLen(_tddist);
    setFltLen(fltlen);
// update electroncs signal time
    updateSignalTime();
// compute initial hit t0 and drift
    updateHitT0(trkt0);
//    std::cout << "creating TrkStrawHit " << this << std::endl;
  }
  
 
  TrkStrawHit::~TrkStrawHit(){
// delete the hit
//    std::cout << "deleting hit " << _theHit << std::endl;
    delete _hittraj;
// ugly trick to keep the base class from trying to delete _TrkDummyHit
    _parentRep=0;  
//    std::cout << "deleted TrkStrawHit " << this << std::endl;
  }

  double
  TrkStrawHit::time() const {
    return strawHit().time();
  }
  
  void
  TrkStrawHit::updateDrift() {
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
// deal with ambiguity updating.  This is a DEPRECATED OPTION, use external ambiguity resolution algorithms instead!!!
    if(_ambigupdate) {
      int iamb = poca().doca() > 0 ? 1 : -1;
      setAmbig(iamb);
    }
// compute the drift time
    double tdrift = strawHit().time() - _hitt0._t0 - _stime;
// find the track direction at this hit
    Hep3Vector tdir = getParentRep()->traj().direction(fltLen());
// convert time to distance.  This computes the intrinsic drift radius error as well
    tcal->TimeToDistance(straw().index(),tdrift,tdir,_t2d);
// Propogate error in t0, using local drift velocity
    double rt0err = _hitt0._t0err*_t2d._vdrift;
    // total hit error is the sum of all
    _toterr = sqrt(_t2d._rdrifterr*_t2d._rdrifterr + rt0err*rt0err + _exterr*_exterr + _penerr*_penerr);
// If the hit is wildly away from the track , disable it
    double rstraw = _straw.getRadius(); 
    if(!physicalDrift(_maxdriftpull)){
      setActivity(false);
      setFlag(driftFail);
    } else {
// otherwise restrict to a physical range
      if (_t2d._rdrift < 0.0){
	_t2d._rdrift = 0.0;
      } else if( _t2d._rdrift > rstraw){
	_t2d._rdrift = rstraw;
      }
    }
  }

  bool
  TrkStrawHit::physicalDrift(double maxchi) const {
    return _t2d._rdrift < _straw.getRadius() + maxchi*_toterr &&
      _t2d._rdrift > -maxchi*_toterr;
  }
  
  void
  TrkStrawHit::updateSignalTime() {
// compute the electronics propagation time.  The convention is that the hit time is measured at the
// FAR END of the wire, as signed by the wire direction.
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    double vwire = tcal->SignalVelocity(straw().index());
    if( poca().status().success()){
      _stime = (straw().getHalfLength()-hitLen())/vwire;
    } else {
// if we're missing poca information, use time division instead
      _stime = (straw().getHalfLength()-_tddist)/vwire;
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
   if( poca().status().success()) {
      status = poca().status();
// update the signal propagation time along the wire
      updateSignalTime();
// update the drift distance using this traj direction
      updateDrift();
// sign drift distance by ambiguity.  Note that an ambiguity of 0 means to ignore the drift
      double residual = poca().doca() - _t2d._rdrift*_iamb;
      setHitResid(residual);
      setHitRms(_toterr);
    } else {
//      cout << "TrkStrawHit:: updateMeasurement() failed" << endl;
      setFlag(updateFail);
      setHitResid(999999);
      setHitRms(999999);
      setActivity(false);
    }
    return status;
  }
  
  void
  TrkStrawHit::hitPosition(Hep3Vector& hpos) const{
    if( poca().status().success() && _iamb!=0){
      Hep3Vector pdir = (trkTraj()->position(fltLen()) - hitTraj()->position(hitLen())).unit();
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
    double inRadius = radius-thick;
    if(pdist>=inRadius)
      pdist = 0.96*inRadius;
    double wallpath =  (sqrt( (radius+pdist)*(radius-pdist) ) -
      sqrt( (inRadius+pdist)*(inRadius-pdist) ));
  // scale for the other dimension
    double cost = tdir.dot(_straw.getDirection());
    if(fabs(cost)<0.999)
      wallpath /= sqrt( (1.0-cost)*(1.0+cost) );
    else
      wallpath = radius;
// use half-length as maximum length
    wallpath = std::min(wallpath,radius);
// test for NAN	
    if(wallpath != wallpath){
      std::cout << "NAN wall" << std::endl;
      wallpath = thick;
   }
    return wallpath;
  }
  
  // compute the pathlength through half the gas , given the drift distance and straw geometry
  double
  TrkStrawHit::gasPath(double pdist,Hep3Vector const& tdir) const {
    double thick = straw().getThickness();
    double radius = straw().getRadius();
    radius-=thick;
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
    if(gaspath != gaspath){
      std::cout << "NAN gas" << endl;
      gaspath = 2*radius;
    }
    return gaspath;
  }

  void TrkStrawHit::print(std::ostream& o) const {
    o<<"------------------- TrkStrawHit -------------------"<<std::endl;
    o<<"istraw "<<_istraw<<std::endl;
    o<<"is active "<<isActive()<<std::endl;
    o<<"hitRms "<<hitRms()<<" weight "<<weight()<<" fltLen "<<fltLen()<<" hitLen "<<hitLen()<<std::endl;
    _strawhit.print(o,true);
    o<<"driftRadius "<<driftRadius()<<" driftRadiusErr "<<driftRadiusErr();
    o<<" hitT0 "<<_hitt0.t0()<<" hitT0err "<<_hitt0.t0Err()<<" t0err "<<t0Err()<<std::endl;
    o<<"ambig "<<_iamb<<" ambig upd "<<_ambigupdate<<std::endl;
    Hep3Vector hpos; hitPosition(hpos);
    o<<"hitPosition "<<hpos<<std::endl;
    o<<"---------------------------------------------------"<<std::endl;
  }

// utility function: this lives in namespace mu2e
  void
  convert(TrkHitVector const& thv, TrkStrawHitVector& tshv) {
    tshv.clear();
    tshv.reserve(thv.size());
    for(auto ith=thv.begin(); ith!=thv.end(); ++ith){
      tshv.push_back(static_cast<TrkStrawHit*>(*ith));
    }
  }

} 
