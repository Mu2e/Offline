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

#include "cetlib_except/coded_exception.h"
#include <algorithm>

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e
{
  TrkStrawHit::TrkStrawHit(StrawResponse::cptr_t strawResponse,
			   const ComboHit& strawhit  , const Straw& straw    , StrawHitIndex index,
			   const TrkT0&    hitt0     , double       fltlen   , double maxdriftpull,
			   double          timeWeight) :
    _strawResponse(strawResponse),
    _combohit(strawhit),
    _straw(straw),
    _index(index),
    _penerr(0.0),
    _toterr(0.0),
    _iamb(0),
    _maxdriftpull(maxdriftpull)
  {
// make sure this ComboHit represents only a single straw hit
    if(_combohit.nStrawHits() != 1 || _combohit.driftEnd() == StrawEnd::unknown)
      throw cet::exception("RECO")<<"mu2e::TrkStrawHit: ComboHit > 1 StrawHit"<< endl;
    // The StrawResponse should be passsed in from outside FIXME!
    Hep3Vector const& wiredir = straw.getDirection();
    Hep3Vector const& mid = straw.getMidPoint();
    // cache the propagation velocity: this depends just on the pulseheight
    _vprop = 2.0*_strawResponse->halfPropV(_combohit.strawId(),1000.0*_combohit.energyDep()); // edep in KeV, FIXME!
    // initialize wire position using time difference
    _wpos = mid +_combohit.wireDist()*wiredir;
// the hit trajectory is defined as a line segment directed along the wire direction starting from the wire center
// ugly conversion to HepPoint FIXME!
    _hittraj = new TrkLineTraj(HepPoint(mid.x(),mid.y(),mid.z()),wiredir,
      timeDiffDist()-timeDiffDistErr(),
      timeDiffDist()+timeDiffDistErr());
    setHitLen(timeDiffDist());
    setFltLen(fltlen);
// update electroncs signal propagation time
    updateSignalTime();
// compute initial hit t0 and drift
//    updateHitT0(hitt0);
    setHitT0(hitt0);
    setHitRms(1.e-6);   // to make sure that the print routine bomb if called from SeedFit
    setActivity(true);
    sett0Weight(timeWeight);
    setTemperature(0.0); // initially no temperature
  }


  TrkStrawHit::~TrkStrawHit(){
// delete the trajectory
    delete _hittraj;
    _parentRep=0;
  }

  double
  TrkStrawHit::driftTime() const {
    return comboHit().time() - _stime - hitT0()._t0;
  }
  

  // bool TrkStrawHit::time( TrkT0& t0 ) {
  //   HitT0 st0;
  //   if (signalPropagationTime(st0)){
  //     t0._t0    = comboHit().time() - st0._t0;
  //     t0._t0err = st0._t0err;
  //     return true;
  //   }else {
  //     return false;
  //   }
  // }


  bool TrkStrawHit::signalPropagationTime( TrkT0& t0 ){
// propagation includes drift and signal propagation along the wire.  First, compute the drift
// time from the distance of closest approach
// correct this for the most recent fit, excluding the effect of this hit on the fit
    double resid, residerr;
    bool retval = this->resid(resid,residerr,true);
    if(retval ) {
      double doca;
    // 0-ambig hit residual is WRT the wire
      if(_iamb == 0)
	doca = fabs(resid);
      else
	doca = _rdrift + _iamb*resid;
    // restrict the range, symmetrically to avoid bias
      double rad       = _straw.getRadius();
      double mint0doca = _strawResponse->Mint0doca(); 
      if(doca > mint0doca && doca < rad-mint0doca){
	// compute phi WRT BField for lorentz drift.
	CLHEP::Hep3Vector trjDir(parentRep()->traj().direction(fltLen()));
	Hep3Vector tperp = trjDir - trjDir.dot(straw().getDirection())*straw().getDirection();
	double phi = tperp.theta();   // This assumes B along z, FIXME!
	// translate the DOCA into a time
	double tdrift = _strawResponse->driftDistanceToTime(_combohit.strawId(), doca, phi);
	double vdrift = _strawResponse->driftInstantSpeed(_combohit.strawId(),doca, phi);
	t0._t0 = tdrift + _stime;
	t0._t0err = residerr/vdrift;// instantaneous velocity to translate the error on the residual
      } else {
	retval = false;
      }
    }
    return retval;
  }

  void
  TrkStrawHit::trackT0Time(double &Htime, double T0flt, const TrkDifPieceTraj* Ptraj, double Vflt){
    throw cet::exception("RECO")<<"mu2e::TrkStrawHit: obsolete function"<< endl;
  }

  void
  TrkStrawHit::updateDrift() {
// compute the drift time
    double tdrift = driftTime();
// find the track direction at this hit
    Hep3Vector tdir = getParentRep()->traj().direction(fltLen());
// convert time to distance.  This computes the intrinsic drift radius error as well

   Hep3Vector tperp = tdir - tdir.dot(straw().getDirection())*straw().getDirection();
   _phi = tperp.theta();
   _rdrift = _strawResponse->driftTimeToDistance(_combohit.strawId(),tdrift,_phi);
   _vdriftinst = _strawResponse->driftInstantSpeed(_combohit.strawId(),fabs(poca().doca()),_phi);
   double vdriftconst = _strawResponse->driftConstantSpeed();
   _rdrifterr = _strawResponse->driftDistanceError(_combohit.strawId(),_rdrift,_phi,fabs(poca().doca()));

// Propogate error in t0, using local drift velocity
    double rt0err = hitT0()._t0err*_vdriftinst;
    // annealing error depends on the 'temperature'
    double exterr = vdriftconst*temperature();
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
    if( poca().status().success()){
      switch (_combohit.driftEnd()) {
	case StrawEnd::cal:
	  _stime = (straw().halfLength()+hitLen())/_vprop;
	  break;
	case StrawEnd::hv:
	  _stime = (straw().halfLength()-hitLen())/_vprop;
	  break;
      }
    } else {
// if we're missing poca information, use time division instead
      switch (_combohit.driftEnd()) {
	case StrawEnd::cal:
	  _stime = (straw().halfLength()+timeDiffDist())/_vprop;
	  break;
	case StrawEnd::hv:
	  _stime = (straw().halfLength()-timeDiffDist())/_vprop;
	  break;
      }
    }
  }

  void
  TrkStrawHit::setAmbig(int newambig){
// if the state changes and the hit is active, warn the rep
    if(isActive() && newambig != _iamb && parentRep() != 0){
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
    _combohit.print(o,true);
    o<<"driftRadius "<<driftRadius()<<" driftRadiusErr "<<driftRadiusErr();
    o<<" hitT0 "<<hitT0().t0()<<" hitT0err "<<hitT0().t0Err()<<" t0err "<<t0Err()<<std::endl;
    o<<"ambig "<<_iamb<< std::endl;
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
