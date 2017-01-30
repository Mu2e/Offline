#include "BTrkData/inc/TrkCaloHit.hh"
// BTrk
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/TrkDifTraj.hh"
#include "BTrk/TrkBase/TrkRep.hh"
#include "BTrk/TrkBase/TrkHit.hh"
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
  TrkCaloHit::TrkCaloHit(const CaloCluster& caloCluster, Hep3Vector &caloClusterPos, 
			 double crystalHalfLength, Hep3Vector const& clusterAxis,
			 const HitT0& hitt0,double fltlen) :
    _caloCluster(caloCluster)
  {

    caloClusterPos.setZ(caloClusterPos.z() + crystalHalfLength);

// the hit trajectory is defined as a line segment directed along the wire direction starting from the wire center
    _hittraj = new TrkLineTraj(HepPoint(caloClusterPos.x(), caloClusterPos.y(), caloClusterPos.z()),
			       clusterAxis, -crystalHalfLength, crystalHalfLength);
    setHitLen(crystalHalfLength);
    setFltLen(fltlen);
// compute initial hit t0 
    updateHitT0(hitt0);
    setActivity(true);
//    std::cout << "creating TrkCaloHit " << this << std::endl;
  }


  TrkCaloHit::~TrkCaloHit(){
// delete the trajectory
    delete _hittraj;
    _parentRep=0;
  }

  double
  TrkCaloHit::time() const {
    return caloCluster().time();
  }

  TrkErrCode
    TrkCaloHit::updateMeasurement(const TrkDifTraj* traj) {
      TrkErrCode status(TrkErrCode::fail);
// find POCA to the wire
    updatePoca(traj);
   if( poca().status().success()) {
      status = poca().status();
// sign drift distance by ambiguity.  Note that an ambiguity of 0 means to ignore the drift
      double residual = poca().doca();//FIXME - _t2d._rdrift*_iamb;
      setHitResid(residual);
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
  TrkCaloHit::hitPosition(Hep3Vector& hpos) const{
    hitTraj()->position(hitLen());
    //hpos = _caloClusterPos;
  }


  void TrkCaloHit::print(std::ostream& o) const {
    o<<"------------------- TrkCaloHit -------------------"<<std::endl;
    // o<<"istraw "<<_istraw<<std::endl;
    // o<<"is active "<<isActive()<<std::endl;
    o<<"hitRms "<<hitRms()<<" weight "<<weight()<<" fltLen "<<fltLen()<<" hitLen "<<hitLen()<<std::endl;
    o<<" hitT0 "<<_hitt0.t0()<<" hitT0err "<<_hitt0.t0Err()<<std::endl;
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
