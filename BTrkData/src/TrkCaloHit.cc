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
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include <algorithm>

#include "TrkReco/inc/TrkUtilities.hh"

using namespace std;
using CLHEP::Hep3Vector;

namespace mu2e
{
  TrkCaloHit::TrkCaloHit(const CaloCluster& caloCluster, Hep3Vector &caloClusterPos, 
			 double crystalHalfLength, Hep3Vector const& clusterAxis,
			 const HitT0& hitt0,double fltlen, double timeWeight, double dtoffset) :
    _caloCluster(caloCluster),
    _dtoffset(dtoffset)
  {

    caloClusterPos.setZ(caloClusterPos.z() + crystalHalfLength);

// the hit trajectory is defined as a line segment directed along the wire direction starting from the wire center
    _hittraj = new TrkLineTraj(HepPoint(caloClusterPos.x(), caloClusterPos.y(), caloClusterPos.z()),
			       clusterAxis, -crystalHalfLength, crystalHalfLength);
    setHitLen(crystalHalfLength);
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

  double
  TrkCaloHit::time() const{
    return caloCluster().time();
  }

  TrkErrCode
  TrkCaloHit::updateMeasurement(const TrkDifTraj* traj) {
    TrkErrCode status(TrkErrCode::fail);
    // find POCA to the wire
    updatePoca(traj);
    if( poca().status().success()) {
      status = poca().status();
      double residual = poca().doca();
      setHitResid(residual);
      double     extErr  = temperature();
      double     totErr  = sqrt(_hitErr*_hitErr + extErr*extErr);
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


  bool 
  TrkCaloHit::signalPropagationTime(double &propTime, double&Doca,
				    double resid    , double &residErr,
				    CLHEP::Hep3Vector trajDirection){
    
    propTime = 0;//FIX ME! 
    residErr = 0.5;//FIX ME!
    return true;
  }
  
  void
  TrkCaloHit::trackT0Time(double& htime, double t0flt, const TrkDifPieceTraj* ptraj, double vflt){
    // compute the flightlength to this hit from z=0 
    CLHEP::Hep3Vector hpos;
    hitPosition(hpos);
    double hflt  = ptraj->zFlight(hpos.z()) - t0flt;
    htime = time() + _dtoffset - hflt/vflt;
  }

  double
  TrkCaloHit::physicalTime() const {
    Hep3Vector trjPos(0);
    hitPosition(trjPos);
    double     CsI_refractiveIndex(1.8);//FIXME! that should come from the geometryhandle
    double     extErr  = temperature();
    double     totErr  = sqrt(_hitErr*_hitErr + extErr*extErr);
    double     caloClX = _caloCluster.cog3Vector().x();
    double     caloClY = _caloCluster.cog3Vector().y();
    double     residx2 = (trjPos.x() - caloClX)*(trjPos.x() - caloClX); 
    double     residy2 = (trjPos.y() - caloClY)*(trjPos.y() - caloClY); 
    double     resid   = sqrt( residx2 + residy2)/(CLHEP::c_light/CsI_refractiveIndex)/totErr;
    
    return resid;
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
