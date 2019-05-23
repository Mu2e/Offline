#include "BTrkData/inc/TrkCaloHit.hh"
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
#include <algorithm>

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e
{
  TrkCaloHit::TrkCaloHit(CaloCluster const& caloCluster, Hep3Vector const& caloClusterPos,
			 double crystalLength, Hep3Vector const& clusterAxis,
			 const HitT0& hitt0,double fltlen, double timeWeight,
			 double hiterr, double terr, double dtoffset) :
    _caloCluster(caloCluster), _clen(crystalLength),
    _dtoffset(dtoffset),
    _hitErr(hiterr) , _tErr(terr)
  {
// the hit trajectory is defined as a line segment directed along the cluster axis 
    _hittraj = new TrkLineTraj(HepPoint(caloClusterPos.x(), caloClusterPos.y(), caloClusterPos.z()),
			       clusterAxis, 0.0, crystalLength);
    setHitLen(0.5*crystalLength); // approximpate
    setFltLen(fltlen);
// compute initial hit t0
    setHitT0(hitt0);
    setActivity(true);

    sett0Weight(timeWeight);
//    std::cout << "creating TrkCaloHit " << this << std::endl;
  }


  TrkCaloHit::~TrkCaloHit(){
// delete the trajectory
    delete _hittraj;
    _parentRep=0;
  }


  //2019-05-02 Gianipez: the following function will change meaning in the near future. FIXME!
  double
  TrkCaloHit::time() const{
    return caloCluster().time()  + _dtoffset; // following Pasha's convention
  }

  // bool 
  // TrkCaloHit::time(HitT0& t0) const{
  //   HitT0 st0;
  //   if (signalPropagationTime(st0)){
  //     t0._t0    = caloCluster().time() - st0._t0 -_dtoffset;
  //     t0._t0err = st0._t0err;
  //     return true;
  //   }else {
  //     return false;
  //   }
  // }

  TrkErrCode
  TrkCaloHit::updateMeasurement(const TrkDifTraj* traj) {
    TrkErrCode status(TrkErrCode::fail);
    // find POCA to the wire
    updatePoca(traj);
    if( poca().status().success()) {
// check the cluster distance to make sure we're on the right loop
      if(hitLen() < _hittraj->lowRange() || hitLen() > 1.5*_hittraj->hiRange()){
	double cost = traj->direction(fltLen()).dot(_hittraj->direction(hitLen()));
	double smax =  0.5*_hittraj->hiRange(); // approximate shwowermax
	double dflt = (hitLen()-smax)/cost;
	setFltLen(fltLen() - dflt);
	setHitLen(smax);
	updatePoca(traj);
      }
      status = poca().status();
      double residual = poca().doca();
      setHitResid(residual);
      double     totErr  = _hitErr; // geometric error is unaffected by temperature
      setHitRms(totErr);
    } else {
//      cout << "TrkCaloHit:: updateMeasurement() failed" << endl;
      setFlag(updateFail);
      setHitResid(999999);
      setHitRms(999999);
      setActivity(false);
    }
    return status;
  }

  void
  TrkCaloHit::hitPosition(CLHEP::Hep3Vector& hpos) const{
    hpos.setX(hitTraj()->position(hitLen()).x());
    hpos.setY(hitTraj()->position(hitLen()).y());
    hpos.setZ(hitTraj()->position(hitLen()).z());
    //hpos = _caloClusterPos;
  }


  bool TrkCaloHit::signalPropagationTime(TrkT0& t0) {
  // compute the light propagation time.
  // light propagation velocity should come from configuration FIXME!
    static const double vlprop =200.0; // mm/nsec  Needs better calibration FIXME!!
    double tlight =0.0;
    if(poca().status().success()){
// time for light to get to SIPMs at the back of the crystals, bounded by crystal length
      double clen = _clen-std::min(_clen,std::max(0.0,poca().flt2()));
      tlight = clen/vlprop;
    }
    t0._t0    =  tlight;
    t0._t0err = _tErr; // intrinsic error on time, used in T0 updating. 
                       //Contribution from the uncertainty of the light propagation is below 100 ps

    return true;//FIXME!
  }

// this function isn't used and needs to be removed FIXME!
  void
  TrkCaloHit::trackT0Time(double& htime, double t0flt, const TrkDifPieceTraj* ptraj, double vflt){
    throw cet::exception("RECO")<<"mu2e::TrkCaloHit: obsolete function"<< endl;
  }

  bool
  TrkCaloHit::isPhysical(double maxchi) const {
    return true;//FIXME!
  }

  void TrkCaloHit::print(std::ostream& o) const {
    o<<"------------------- TrkCaloHit -------------------"<<std::endl;
    // o<<"istraw "<<_istraw<<std::endl;
    // o<<"is active "<<isActive()<<std::endl;
    o<<"hitRms "<<hitRms()<<" weight "<<weight()<<" fltLen "<<fltLen()<<" hitLen "<<hitLen()<<std::endl;
    o<<" hitT0 "<<hitT0().t0()<<" hitT0err "<<hitT0().t0Err()<<std::endl;
    Hep3Vector hpos; hitPosition(hpos);
    o<<"hitPosition "<<hpos<<std::endl;
    o<<"---------------------------------------------------"<<std::endl;
  }

// utility function: this lives in namespace mu2e
  void
  convert(TrkHitVector const& thv, TrkCaloHitVector& tshv) {
    tshv.clear();
    tshv.reserve(thv.size());
    for(auto ith=thv.begin(); ith!=thv.end(); ++ith){
      TrkCaloHit* tsh = dynamic_cast<TrkCaloHit*>(*ith);
      if(tsh != 0) tshv.push_back(tsh);
    }
  }

}
