//
// BaBar hit object corresponding to a single cell hit
//
// $Id: TrkCellHit.cc,v 1.2 2012/12/04 00:51:27 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:27 $
//
// Original author G. Tassielli
//

#include "BaBar/BaBar.hh"
#include "KalmanTestsI/inc/TrkCellHit.hh"
//#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/TrkPoca.hh"
#include "TrkBase/TrkDifTraj.hh"
//#include "TrkBase/TrkDetElemId.hh"
#include "TrkBase/TrkRep.hh"
//#include "TrkBase/HelixParams.hh"

//geometry
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrationsI.hh"
//#include "ConditionsService/inc/TrackerCalibrations.hh"

#include "cetlib/pow.h"

#include <iostream>

using namespace std;

namespace mu2e
{
/*
/// Material information, BaBar style
  MatDBInfo* TrkCellHit::matdbinfo(){
    static MatDBInfo mat;
    return &mat;
  }*/
  DetStrawHitType* TrkCellHit::cgtype(std::string gasMat){
    static DetStrawHitType instance(matdbinfo(),gasMat.c_str());
    return &instance;
  }
  DetStrawHitType* TrkCellHit::fwtype(std::string fwMat){
    static DetStrawHitType instance(matdbinfo(),fwMat.c_str());
    return &instance;
  }
  DetStrawHitType* TrkCellHit::swtype(std::string swMat){
    static DetStrawHitType instance(matdbinfo(),swMat.c_str());
    return &instance;
  }

  TrkCellHit::TrkCellHit(const StrawHit& strawhit, const Straw& straw, unsigned istraw,
    const TrkT0& trkt0,double fltlen,double exterr,double maxdriftpull) :
    TrkStrawHit(strawhit, straw, istraw,
        trkt0, fltlen, exterr, maxdriftpull),
        _cgelem(cgtype(),this),
        _fwelem_b(fwtype(),this),
        _fwelem_s(fwtype(),this,1),
        _fwelem_t(fwtype(),this,2),
        _swelem(swtype(),this),
        _evalCellPath(true),
        _entryDltPath(0.0),
        _exitDltPath(0.0),
        _growingRad(false),
        _hitOutOfCut(false),
        _hitSWire(false),
        _hitFWire_t(false),
        _hitFWire_s(false),
        _hitFWire_b(false),
        _pathlInSWire(0.0),
        _pathlInFWire_t(0.0),
        _pathlInFWire_s(0.0),
        _pathlInFWire_b(0.0)
  {
          const Tracker& tracker = getTrackerOrThrow();
          const mu2e::ITracker &itracker = static_cast<const mu2e::ITracker&>( tracker );
          _itwp = itracker.getCellGeometryHandle();
  }

  TrkCellHit::TrkCellHit(const StrawHit& strawhit, const Straw& straw, unsigned istraw,
    const TrkT0& trkt0,double fltlen,double exterr,double maxdriftpull, std::string matOnly) :
    TrkStrawHit(strawhit, straw, istraw,
        trkt0, fltlen, exterr, maxdriftpull),
        _cgelem(cgtype(matOnly),this),
        _fwelem_b(fwtype(matOnly),this),
        _fwelem_s(fwtype(matOnly),this,1),
        _fwelem_t(fwtype(matOnly),this,2),
        _swelem(swtype(matOnly),this),
        _evalCellPath(true),
        _entryDltPath(0.0),
        _exitDltPath(0.0),
        _growingRad(false),
        _hitOutOfCut(false),
        _hitSWire(false),
        _hitFWire_t(false),
        _hitFWire_s(false),
        _hitFWire_b(false),
        _pathlInSWire(0.0),
        _pathlInFWire_t(0.0),
        _pathlInFWire_s(0.0),
        _pathlInFWire_b(0.0)
  {
          const Tracker& tracker = getTrackerOrThrow();
          const mu2e::ITracker &itracker = static_cast<const mu2e::ITracker&>( tracker );
          _itwp = itracker.getCellGeometryHandle();
  }

  TrkCellHit::TrkCellHit(const TrkCellHit& other, TrkRep* rep) :
                  TrkStrawHit(other, rep),
                  _cgelem(cgtype(),this),
                  _fwelem_b(fwtype(),this),
                  _fwelem_s(fwtype(),this,1),
                  _fwelem_t(fwtype(),this,2),
                  _swelem(swtype(),this),
                  _evalCellPath(true),
                  _entryDltPath(0.0),
                  _exitDltPath(0.0),
                  _growingRad(false),
                  _hitOutOfCut(false),
                  _hitSWire(false),
                  _hitFWire_t(false),
                  _hitFWire_s(false),
                  _hitFWire_b(false),
                  _pathlInSWire(0.0),
                  _pathlInFWire_t(0.0),
                  _pathlInFWire_s(0.0),
                  _pathlInFWire_b(0.0)
  {
    //std::cout << "creating TrkCellHit copy " << this << std::endl;
          const Tracker& tracker = getTrackerOrThrow();
          const mu2e::ITracker &itracker = static_cast<const mu2e::ITracker&>( tracker );
          _itwp = itracker.getCellGeometryHandle();

  }

//  TrkCellHit::TrkCellHit(const TrkCellHit& other, const TrkDifTraj* trkTraj) :
//                  TrkStrawHit(other, trkTraj),
//                  _cgelem(cgtype(),this),
//                  _fwelem_b(fwtype(),this),
//                  _fwelem_s(fwtype(),this,1),
//                  _fwelem_t(fwtype(),this,2),
//                  _swelem(swtype(),this),
//                  _evalCellPath(true),
//                  _entryDltPath(0.0),
//                  _exitDltPath(0.0),
//                  _growingRad(false),
//                  _hitOutOfCut(false),
//                  _hitSWire(false),
//                  _hitFWire_t(false),
//                  _hitFWire_s(false),
//                  _hitFWire_b(false),
//                  _pathlInSWire(0.0),
//                  _pathlInFWire_t(0.0),
//                  _pathlInFWire_s(0.0),
//                  _pathlInFWire_b(0.0)
//  {
//    //std::cout << "creating TrkCellHit copy " << this << std::endl;
//          const Tracker& tracker = getTrackerOrThrow();
//          const mu2e::ITracker &itracker = static_cast<const mu2e::ITracker&>( tracker );
//          _itwp = itracker.getCellGeometryHandle();
//
//  }

  TrkCellHit::~TrkCellHit(){

  }

  TrkCellHit*
  TrkCellHit::clone(TrkRep* parentRep, const TrkDifTraj* trkTraj) const {
    return new TrkCellHit(*this, parentRep);
  }

  CellGeometryHandle* TrkCellHit::cellHandle() const {
          _itwp->SelectCellDet(strawHit().strawIndex().asUint());
          return _itwp;
  }

  int TrkCellHit::GetSuperLayer() const {
          _itwp->SelectCellDet(strawHit().strawIndex().asUint());
          return _itwp->GetSuperLayer();
  }
  int TrkCellHit::GetCelRing() const {
          _itwp->SelectCellDet(strawHit().strawIndex().asUint());
          return _itwp->GetCelRing();
  }
  int TrkCellHit::GetCell() const {
          _itwp->SelectCellDet(strawHit().strawIndex().asUint());
         return _itwp->GetWire();
  }

/*
  double
  TrkCellHit::time() const {
    return strawHit().time();
  }
*/
  void
  TrkCellHit::updateDrift() {
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
// deal with ambiguity updating.  This is a DEPRECATED OPTION, use external ambiguity resolution algorithms instead!!!
    if(_ambigupdate) {
      int iamb = _poca->doca() > 0 ? 1 : -1;
      setAmbig(iamb);
    }
// compute the drift time
    double tdrift = strawHit().time() - _hitt0._t0 - _stime;
// find the track direction at this hit
    CLHEP::Hep3Vector tdir = getParentRep()->traj().direction(fltLen());
// convert time to distance.  This computes the intrinsic drift radius error as well
    tcal->TimeToDistance(strawHit().strawIndex(),tdrift,tdir,_t2d);
// Propogate error in t0, using local drift velocity
    double rt0err = _hitt0._t0err*_t2d._vdrift;
    // total hit error is the sum of all
    _toterr = sqrt(_t2d._rdrifterr*_t2d._rdrifterr + rt0err*rt0err + _exterr*_exterr + _penerr*_penerr);
// If the hit is wildly away from the track , disable it
    double rstraw = cellHandle()->GetCellRad();
    if( _t2d._rdrift - rstraw > _maxdriftpull*_toterr ||
      _t2d._rdrift < -_maxdriftpull*_toterr){
        setUsability(10);//-10
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
  TrkCellHit::updateSignalTime() {
// compute the electronics propagation time.  The convention is that the hit time is measured at the
// FAR END of the wire, as signed by the wire direction.
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    double vwire = tcal->SignalVelocity(strawHit().strawIndex());
    _itwp->SelectCellDet(strawHit().strawIndex().asUint());
    if(_poca != 0 && _poca->status().success()){
            if (hitLen()>_itwp->GetCellHalfLength()) {
                    _stime = 0.0;
            } else if (hitLen()<-_itwp->GetCellHalfLength()) {
                    _stime = 2.0*_itwp->GetCellHalfLength()/vwire;
            } else {
                    _stime = (_itwp->GetCellHalfLength()-hitLen())/vwire;
            }
    } else {
// if we're missing poca information, use time division instead
      if (_tddist>_itwp->GetCellHalfLength()) {
              _stime = 0.0;
      } else if (_tddist<-_itwp->GetCellHalfLength()) {
              _stime = 2.0*_itwp->GetCellHalfLength()/vwire;
      } else {
              _stime = (_itwp->GetCellHalfLength()-_tddist)/vwire;
      }
    }
  } 
/*
  void 
  TrkCellHit::setAmbig(int newambig){
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
  TrkCellHit::updateMeasurement(const TrkDifTraj* traj) {
    TrkErrCode status(TrkErrCode::fail);
// find POCA to the wire
    updatePoca(traj);
   if(_poca != 0 && _poca->status().success()) {
      status = _poca->status();
// update the signal propagation time along the wire
      updateSignalTime();
// update the drift distance using this traj direction
      updateDrift();
// sign drift distance by ambiguity.  Note that an ambiguity of 0 means to ignore the drift
      double residual = _poca->doca() - _t2d._rdrift*_iamb;
      setHitResid(residual);
      setHitRms(_toterr);
    } else {
      cout << "TrkCellHit:: updateMeasurement() failed" << endl;
      setHitResid(999999);
      setHitRms(999999);
      setUsability(0);
    }
    return status;
  }
  
  void
  TrkCellHit::hitPosition(CLHEP::Hep3Vector& hpos) const{
    if(_poca != 0 && _poca->status().success() && _iamb!=0){
      CLHEP::Hep3Vector pdir = (trkTraj()->position(fltLen()) - hitTraj()->position(hitLen())).unit();
      hpos = _wpos + pdir*_t2d._rdrift*_iamb;
    } else {
      hpos = _wpos;
    }
  }
*/

  // compute the total pathlength through one Cell, given the drift distance and Cell geometry
const double TrkCellHit::cellPath(double pftl,Hep3Vector const& tdir, HepPoint const& tpos, double omega, double cosDip, double angle){

        //cout<<"pftl "<<pftl<<" cell MaxRad "<<cellHandle()->GetCellRad()<<" insideRad "<<cellHandle()->GetCellInsideRad()<<endl;

        _hitOutOfCut=false;
        _centerDlPath=0.0;
        const Hep3Vector &wDir = cellHandle()->GetWireDirection();
        Hep3Vector tdirPw = tdir - (tdir.dot(wDir))*wDir;
        //cout<<"tdir "<<tdir<<" tdirPw "<<tdirPw<<endl;

        double ttpos[3], tposCellRef[3];
        ttpos[0]=tpos.x();
        ttpos[1]=tpos.y();
        ttpos[2]=tpos.z();
        _itwp->Global2Local(ttpos,tposCellRef);
        CLHEP::Hep3Vector wCenter;
        _itwp->WirePosAtZ(tpos.z(),wCenter);

        //cout<<"cell center at z "<<tpos.z()<<" = "<<wCenter<<endl;
        //cout<<"global pos "<<tpos<<" relative to cell pos("<<tposCellRef[0]<<","<<tposCellRef[1]<<","<<tposCellRef[2]<<")"<<endl;
        //cout<<"wire epsilon "<<_itwp->GetWireEpsilon()<<endl;
        double maxCellDim = _itwp->GetCellInsideRad();
        double cutMaxCellDim = 1.2*maxCellDim;
        double maxDelta = 2.0*maxCellDim;
        double maxCut = maxDelta*1.5;

        double deltaX1(0.0);
        double deltaY1(0.0);
        double deltaX2(0.0);
        double deltaY2(0.0);

        double angleCellDir = wCenter.azimAngle(/*tdir*/tdirPw);
        double angleAdjast = angleCellDir;
        //cout<<"angleCellDir "<<angleCellDir<<endl;
        double xSign = 1.0, ySign=-1.0;
        _growingRad = true;
        if (angleCellDir<0.0) {
                xSign = -1.0;
                angleAdjast = -angleCellDir;
        }
        if (angleCellDir<-CLHEP::halfpi || angleCellDir>CLHEP::halfpi) {
                ySign= 1.0;
                angleAdjast = CLHEP::pi -xSign*angleCellDir;
                _growingRad = false;
        }

        double tmpMaxCell = fabs(tposCellRef[0]);
        if (tmpMaxCell>maxCellDim && tmpMaxCell<cutMaxCellDim) { maxCellDim=tmpMaxCell; }
        tmpMaxCell = fabs(tposCellRef[1]);
        if (tmpMaxCell>maxCellDim && tmpMaxCell<cutMaxCellDim) { maxCellDim=tmpMaxCell; }

        deltaX1 =  -maxCellDim + xSign*tposCellRef[0];
        deltaY1 =  -maxCellDim + ySign*tposCellRef[1];
        deltaX2 =   maxCellDim + xSign*tposCellRef[0];
        deltaY2 =   maxCellDim + ySign*tposCellRef[1];
        if (deltaX1>0.0) { deltaX1=0.0; }
        else if (deltaX1<-maxDelta) { deltaX1=-maxDelta; }
        if (deltaY1>0.0) { deltaY1=0.0; }
        else if (deltaY1<-maxDelta) { deltaY1=-maxDelta; }
        if (deltaX2<0.0) { deltaX2=0.0; }
        else if (deltaX2> maxDelta) { deltaX2= maxDelta; }
        if (deltaY2<0.0) { deltaY2=0.0; }
        else if (deltaY2> maxDelta) { deltaY2= maxDelta; }

        double deltaS1(0.0), deltaS2(0.0);

        //cout<<"deltaX1 "<<deltaX1<<" deltaY1 "<<deltaY1<<" deltaX2 "<<deltaX2<<" deltaY2 "<<deltaY2<<endl;

//        if (omega>1e-4 && omega<0.1 ) {
//                cosDip=fabs(cosDip);
//                cout<<"hit curvature "<<omega<<" cosDip "<<cosDip<<endl;
//                double sToAngle = cosDip*omega;
//                double RcosAngle0 = cos(angle);
//                double RsinAngle0 = sqrt(1.0-RcosAngle0*RcosAngle0)/omega;
//                RcosAngle0/=omega;
//                double deltaStmp(0.0);
//                cout<<"sToAngle "<<sToAngle<<" angle "<<angle<<" RcosAngle0 "<<RcosAngle0<<" RsinAngle0 "<<RsinAngle0<<endl;
//                cout<<"----------------"<<endl;
//                double tmpDelta = RsinAngle0*RsinAngle0/*-*/+2.0*RcosAngle0*deltaX1;
//                if (tmpDelta>0.0) {
//                        tmpDelta = sqrt(tmpDelta);
//                        //deltaS1 = std::min(fabs((-RsinAngle0-tmpDelta)/RcosAngle0),fabs((-RsinAngle0+tmpDelta)/RcosAngle0));
//                        deltaS1 = (-RsinAngle0+tmpDelta)/RcosAngle0;
//                        cout<<"deltaT1 : "<<deltaS1<<endl;
//                        deltaS1 /= sToAngle;
//                        cout<<"deltaS1 : "<<deltaS1<<endl;
//                        cout<<"--- : "<<((-RsinAngle0-tmpDelta)/RcosAngle0)/sToAngle<<" or "
//                                        <<((-RsinAngle0+tmpDelta)/RcosAngle0)/sToAngle<<endl;
//                } else {
//                        deltaS1 = 0.0;//1.e5;
//                }
//
//                tmpDelta = RcosAngle0*RcosAngle0-2.0*RsinAngle0*deltaY1;
//                if (tmpDelta>0.0) {
//                        tmpDelta = sqrt(tmpDelta);
//                        //deltaStmp = std::min(fabs((RcosAngle0-tmpDelta)/RsinAngle0),fabs((RcosAngle0+tmpDelta)/RsinAngle0));
//                        deltaStmp = (RcosAngle0-tmpDelta)/RsinAngle0;
//                        cout<<"deltaTtmp : "<<deltaStmp<<endl;
//                        deltaStmp /= sToAngle;
//                        cout<<"deltaStmp : "<<deltaStmp<<endl;
//                        cout<<"--- : "<<((RcosAngle0-tmpDelta)/RsinAngle0)/sToAngle<<" or "
//                                        <<((RcosAngle0+tmpDelta)/RsinAngle0)/sToAngle<<endl;
//                } else {
//                        deltaStmp = 0.0;//1.e5;
//                }
//                if (normalDire) {
//                        deltaS1 = std::max(deltaS1,deltaStmp);
//                } else {
//                        deltaS1 = std::min(deltaS1,deltaStmp);
//                }
//                //if (deltaS1>0.999e5) deltaS1=0.0;
//
//                tmpDelta = RsinAngle0*RsinAngle0/*-*/+2.0*RcosAngle0*deltaX2;
//                if (tmpDelta>0.0) {
//                        tmpDelta = sqrt(tmpDelta);
//                        //deltaS2 = std::min(fabs((-RsinAngle0-tmpDelta)/RcosAngle0),fabs((-RsinAngle0+tmpDelta)/RcosAngle0));
//                        deltaS2 = (-RsinAngle0+tmpDelta)/RcosAngle0;
//                        cout<<"deltaT2 : "<<deltaS2<<endl;
//                        deltaS2 /= sToAngle;
//                        cout<<"deltaS2 : "<<deltaS2<<endl;
//                        cout<<"--- : "<<((-RsinAngle0-tmpDelta)/RcosAngle0)/sToAngle<<" or "
//                                        <<((-RsinAngle0+tmpDelta)/RcosAngle0)/sToAngle<<endl;
//                } else {
//                        deltaS2 = 0.0;//1.e5;
//                }
//
//                tmpDelta = RcosAngle0*RcosAngle0-2.0*RsinAngle0*deltaY2;
//                if (tmpDelta>0.0) {
//                        tmpDelta = sqrt(tmpDelta);
//                        //deltaStmp = std::min(fabs((RcosAngle0-tmpDelta)/RsinAngle0),fabs((RcosAngle0+tmpDelta)/RsinAngle0));
//                        deltaStmp = (RcosAngle0-tmpDelta)/RsinAngle0;
//                        cout<<"deltaTtmp : "<<deltaStmp<<endl;
//                        deltaStmp /= sToAngle;
//                        cout<<"deltaStmp : "<<deltaStmp<<endl;
//                        cout<<"--- : "<<((RcosAngle0-tmpDelta)/RsinAngle0)/sToAngle<<" or "
//                                        <<((RcosAngle0+tmpDelta)/RsinAngle0)/sToAngle<<endl;
//                } else {
//                        deltaStmp = 0.0;//1.e5;
//                }
//
//                if (normalDire) {
//                        deltaS2 = std::min(deltaS2,deltaStmp);
//                } else {
//                        deltaS2 = std::max(deltaS2,deltaStmp);
//                }
//                //if (deltaS2>0.999e5) deltaS2=0.0;
//                cout<<"----------------"<<endl;
//
//        } else {
//                cout<<"linear approximation"<<endl;

                //double invCosTh = sqrt( cet::sum_of_squares(tdirPw.x(), tdirPw.y()) );
                //double invSinTh = invCosTh/fabs(tdirPw.y()); //the wire frame is 90 deg rotated
                //invCosTh /= fabs(tdirPw.x());
                //double invCosTh

        {
                //cout<<"angleAdjast "<<angleAdjast<<endl;
                double invCosTh = cos(angleAdjast);
                double invSinTh = sqrt(1.0-invCosTh*invCosTh);
                invCosTh = 1.0/invCosTh;
                invSinTh = 1.0/invSinTh;
                //cout<<"invCosTh "<<invCosTh<<" invSinTh "<<invSinTh<<endl;

                deltaX1*=invSinTh; deltaY1*=invCosTh;
                deltaX2*=invSinTh; deltaY2*=invCosTh;

                deltaS1 = std::max(deltaX1, deltaY1);
                if (deltaS1>0) {deltaS1=0.0;}
                else if (deltaS1<-maxDelta) {deltaS1=-maxDelta;}
                deltaS2 = std::min(deltaX2, deltaY2);
                if (deltaS2<0) {deltaS2=0.0;}
                else if (deltaS2>maxDelta) {deltaS2=maxDelta;}
                _centerDlPath = -tposCellRef[1]*invCosTh;

                //cout<<"X1 "<<tposCellRef[0]-xSign*deltaS1/invSinTh<<" Y1 "<<tposCellRef[1]-ySign*deltaS1/invCosTh<<endl;
                //cout<<"X2 "<<tposCellRef[0]-xSign*deltaS2/invSinTh<<" Y2 "<<tposCellRef[1]-ySign*deltaS2/invCosTh<<endl;

                //cout<<"deltaS1 "<<deltaS1<<" deltaS2 "<<deltaS2<<endl;
                if ((deltaS2-deltaS1)>maxCut) {
                        deltaS1 = -0.5*maxCut;
                        deltaS2 =  0.5*maxCut;
                        _centerDlPath = 0.0;
                        _hitOutOfCut=true;
                }

                double cost = tdir.dot(wDir);
                double invsint(0.0);
                if(fabs(cost)<0.999) {
                        invsint = 1.0/sqrt( (1.0-cost)*(1.0+cost) );
                        deltaS1 *= invsint;
                        deltaS2 *= invsint;
                        _centerDlPath *= invsint;
                }
        }
        //cout<<"pos "<<tpos<<" dir "<<tdir<<" tdir.phi() "<<tdir.phi()/*" tdirPw "<<tdirPw<<" tdirPw.phi() "<<tdirPw.phi()*/<<endl;
        //cout<<"deltaS1 "<<deltaS1<<" deltaS2 "<<deltaS2<<endl;

        _entryDltPath =  deltaS1;
        _exitDltPath  =  deltaS2;

        //cout<<"entryPath "<<_entryDltPath<<" exitPath "<<_exitDltPath<<endl;
        _evalCellPath=false;

    //    else
    //      wallpath = radius;
    //// use half-length as maximum length
    //    wallpath = std::min(wallpath,radius);
    //// test for NAN
    ////    if(wallpath != wallpath){
    /////      std::cout << "non wall" << std::endl;
    ////   }
        return (_exitDltPath-_entryDltPath);
}

// compute the pathlength through one wall of the straw, given the drift distance and straw geometry
  const double
  TrkCellHit::fieldWirePath(double pftl,Hep3Vector const& tdir, HepPoint const& tpos, double omega, double cosDip, double angle) {

          if(_evalCellPath) cellPath(pftl, tdir, tpos, omega, cosDip, angle);
          return 0.0;
  }

  // compute the pathlength through one wall of the straw, given the drift distance and straw geometry
    const double
    TrkCellHit::senseWirePath(double pftl,Hep3Vector const& tdir, HepPoint const& tpos, double omega, double cosDip, double angle) {

            if(_evalCellPath) cellPath(pftl, tdir, tpos, omega, cosDip, angle);
            return _centerDlPath;
    }

  // compute the pathlength through half the gas , given the drift distance and straw geometry
  const double
  TrkCellHit::cellGasPath(double pftl,Hep3Vector const& tdir, HepPoint const& tpos, double omega, double cosDip, double angle) {
          double gaspath = 0.0;
          if(_evalCellPath) { gaspath=cellPath(pftl, tdir, tpos, omega, cosDip, angle); }
          else { gaspath =(_exitDltPath-_entryDltPath); }
          _evalCellPath=true;
          return gaspath;
  }

}

