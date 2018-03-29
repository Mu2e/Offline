//
// BaBar hit object corresponding to a single straw hit
//
// Original author David Brown, LBNL
//
#include "BTrkData/inc/TrkStrawHit.hh"
// BTrk
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/TrkDifTraj.hh"
#include "BTrk/TrkBase/TrkDifPieceTraj.hh"
#include "BTrk/TrkBase/TrkRep.hh"
#include "BTrk/TrkBase/TrkHit.hh"
//
#include "CLHEP/Vector/ThreeVector.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "TrackerConditions/inc/StrawResponse.hh"

#include "TrkReco/inc/TrkUtilities.hh"

#include <algorithm>

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e
{
  TrkStrawHit::TrkStrawHit(const StrawHit& strawhit  , const Straw& straw    , StrawHitIndex index,
			   const TrkT0&    hitt0     , double       fltlen   , double   exterr, double maxdriftpull, 
			   double          timeWeight, double       minT0doca) :
    _strawhit(strawhit),
    _straw(straw),
    _index(index),
    _penerr(0.0),
    _toterr(0.0),
    _iamb(0),
    _enduse(earliest), // default should be set by a configuration object FIXME!
    _ambigupdate(false),
    _maxdriftpull(maxdriftpull),
    _mint0doca(minT0doca)
  {
// The position information should come from the StrawHitPosition collection, FIXME!!! 
    ConditionsHandle<StrawResponse> srep = ConditionsHandle<StrawResponse>("ignored");
    float dw, dwerr;
    srep->wireDistance(strawhit,straw.getHalfLength(),dw,dwerr);
    _wpos = straw.getMidPoint()+dw*straw.getDirection();
    _tddist = dw; 
    _tddist_err = dwerr;
    Hep3Vector const& wiredir = straw.getDirection();
    Hep3Vector const& mid = straw.getMidPoint();
// the hit trajectory is defined as a line segment directed along the wire direction starting from the wire center
    _hittraj = new TrkLineTraj(HepPoint(mid.x(),mid.y(),mid.z()),wiredir,_tddist-_tddist_err,_tddist+_tddist_err);
    setHitLen(_tddist);
    setFltLen(fltlen);
// update electroncs signal time
    updateSignalTime();
// compute initial hit t0 and drift
//    updateHitT0(hitt0);
    setHitT0(hitt0);
    setActivity(true);
//    std::cout << "creating TrkStrawHit " << this << std::endl;
    sett0Weight(timeWeight);
    setTemperature(exterr);
  }


  TrkStrawHit::~TrkStrawHit(){
// delete the trajectory
    delete _hittraj;
    _parentRep=0;
  }

  double
  TrkStrawHit::driftTime(StrawEnd end) const {
    return strawHit().time(end) - hitT0()._t0 - _stime[end];
  }

  double
  TrkStrawHit::driftTime() const {
    double tdrift;
    if(_enduse==earliest){
      tdrift = std::min(strawHit().time(TrkTypes::cal) - _stime[TrkTypes::cal],
	  strawHit().time(TrkTypes::hv) - _stime[TrkTypes::hv] )
	- hitT0()._t0;
    } else if(_enduse ==both){
      tdrift = 0.5*(strawHit().time(TrkTypes::cal) - _stime[TrkTypes::cal] +
	  strawHit().time(TrkTypes::hv) - _stime[TrkTypes::hv] )
	- hitT0()._t0;
    } else {
      tdrift = strawHit().time(static_cast<TrkTypes::End>(_enduse)) - _stime[_enduse] - hitT0()._t0;
    }
    return tdrift;
  }

  bool
  TrkStrawHit::signalPropagationTime(double &PropTime, double &Doca      , 
				     double Resid    , double &ResidError, 
				     CLHEP::Hep3Vector TrajDirection){
    // update the drift distance using this traj direction
    //    updateDrift();

    ConditionsHandle<StrawResponse> srep = ConditionsHandle<StrawResponse>("ignored");
    // convert this to a distance to the wire
    double driftRadius = _rdrift;
    Doca = (Resid + driftRadius*_iamb);
    if(_iamb == 0)
      Doca = fabs(Doca);
    else
      Doca *= _iamb;
    // restrict the range, symmetrically to avoid bias
    double rad       = _straw.getRadius();
    double mint0doca = srep->Mint0doca(); 
    if(Doca > mint0doca && Doca < rad-mint0doca){
      // translate the DOCA into a time
      Hep3Vector tperp = TrajDirection - TrajDirection.dot(straw().getDirection())*straw().getDirection();
      double phi = tperp.theta(); 
      double tdrift = srep->driftDistanceToTime(straw().index(), Doca, phi);
      double vdrift = srep->driftInstantSpeed(straw().index(),Doca, phi);
      PropTime    = tdrift;
      switch(_enduse) {
	case cal: case hv:
	  PropTime += _stime[_enduse];
	  break;
	case earliest:
	  PropTime += std::min(_stime[0],_stime[1]);
	  break;
	case both :
	  PropTime += 0.5*(_stime[0] + _stime[1]);
	  break;
      }
      ResidError /= vdrift;
      return true;
    } else {
      PropTime = 0;
      return false;
    }
  }
  
  void 
  TrkStrawHit::trackT0Time(double &Htime, double T0flt, const TrkDifPieceTraj* Ptraj, double Vflt){
    ConditionsHandle<StrawResponse> srep = ConditionsHandle<StrawResponse>("ignored");
    // compute the flightlength to this hit from z=0 (can be negative)
    double pz   = _wpos.z();
    double hflt = Ptraj->zFlight(pz) - T0flt;
    // Use this to estimate the time for the track to reaches this hit from z=0
    double tprop = hflt/Vflt;
    // estimate signal propagation time on the wire assuming the middle (average)
    double vwire = srep->halfPropV(straw().index(),strawHit().energyDep()*1000.)*2; //FIXME
    double teprop = _straw.getHalfLength()/vwire;
    // correct the measured time for these effects: this gives the aveage time the particle passed this straw, WRT
    // when the track crossed Z=0
    // assume the average drift time is half the maximum drift distance.  This is a poor approximation, but good enough for now
    // for crude estimates, we only need 1 d2t function
    CLHEP::Hep3Vector zdir(0.0,0.0,1.0);
    double phi = 0; // FIXME should default phi be 0?
    double tdrift = srep->driftDistanceToTime(straw().index(), 0.5*straw().getRadius(), phi);
    Htime = _strawhit.time() - tprop - teprop - tdrift;
  }


  void
  TrkStrawHit::updateDrift() {
    ConditionsHandle<StrawResponse> srep = ConditionsHandle<StrawResponse>("ignored");
// deal with ambiguity updating.  This is a DEPRECATED OPTION, use external ambiguity resolution algorithms instead!!!
    if(_ambigupdate) {
      int iamb = poca().doca() > 0 ? 1 : -1;
      setAmbig(iamb);
    }
// compute the drift time
    double tdrift = driftTime();
// find the track direction at this hit
    Hep3Vector tdir = getParentRep()->traj().direction(fltLen());
// convert time to distance.  This computes the intrinsic drift radius error as well

   Hep3Vector tperp = tdir - tdir.dot(straw().getDirection())*straw().getDirection();
   double phi = tperp.theta(); 
   _rdrift = srep->driftTimeToDistance(straw().index(),tdrift,phi);
   _vdriftinst = srep->driftInstantSpeed(straw().index(),_rdrift,phi);
   _rdrifterr = srep->driftDistanceError(straw().index(),_rdrift,phi,fabs(poca().doca()));

// Propogate error in t0, using local drift velocity
    double rt0err = hitT0()._t0err*_vdriftinst;
    // annealing error depends on the 'temperature'
    double exterr = _vdriftinst*temperature();
    // total hit error is the sum of all
    _toterr = sqrt(_rdrifterr*_rdrifterr + rt0err*rt0err + exterr*exterr + _penerr*_penerr);
// If the hit is wildly away from the track , disable it
    double rstraw = _straw.getRadius();
    if(!isPhysical(_maxdriftpull)){
      setActivity(false);
      setFlag(driftFail);
    } else {
// otherwise restrict to a physical range
      if (_rdrift < 0.0){
        _rdrift = 0.0;
      } else if( _rdrift > rstraw){
        _rdrift = rstraw;
      }
    }
  }

  bool TrkStrawHit::isPhysical(double maxchi) const {
    return _rdrift < _straw.getRadius() + maxchi*_toterr &&
      _rdrift > -maxchi*_toterr;
  }

  void
  TrkStrawHit::updateSignalTime() {
// compute the electronics propagation time for the 2 ends.
// note: the wire direction points from cal to HV
    ConditionsHandle<StrawResponse> srep = ConditionsHandle<StrawResponse>("ignored");
    double vwire = srep->halfPropV(straw().index(),strawHit().energyDep()*1000.)*2; //FIXME
    if( poca().status().success()){
      _stime[TrkTypes::cal] = (straw().getHalfLength()+hitLen())/vwire;
      _stime[TrkTypes::hv] = (straw().getHalfLength()-hitLen())/vwire;
    } else {
// if we're missing poca information, use time division instead
      _stime[TrkTypes::cal] = (straw().getHalfLength()+_tddist)/vwire;
      _stime[TrkTypes::hv] = (straw().getHalfLength()-_tddist)/vwire;
    }
  }

  void
  TrkStrawHit::setAmbig(int newambig){
// if the state changes and the hit is active, warn the rep
    if(isActive() && newambig != _iamb){
      parentRep()->setCurrent(false);
//      std::cout << "changing hit ambiguity " << std::endl;
    }
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
      double residual = poca().doca() - _rdrift*_iamb;
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
      hpos = _wpos + pdir*_rdrift*_iamb;
    } else {
      hpos = _wpos;
    }
  }


  void TrkStrawHit::print(std::ostream& o) const {
    o<<"------------------- TrkStrawHit -------------------"<<std::endl;
    o<<"straw hit "<<_index<<std::endl;
    o<<"is active "<<isActive()<<std::endl;
    o<<"hitRms "<<hitRms()<<" weight "<<weight()<<" fltLen "<<fltLen()<<" hitLen "<<hitLen()<<std::endl;
    _strawhit.print(o,true);
    o<<"driftRadius "<<driftRadius()<<" driftRadiusErr "<<driftRadiusErr();
    o<<" hitT0 "<<hitT0().t0()<<" hitT0err "<<hitT0().t0Err()<<" t0err "<<t0Err()<<std::endl;
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
      TrkStrawHit* tsh = dynamic_cast<TrkStrawHit*>(*ith);
      if(tsh != 0) tshv.push_back(tsh);
    }
  }

}
