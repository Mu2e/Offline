//--------------------------------------------------------------------------
// Name:
//   DetStrawGasElem: Dummy class to represent a straw element.
//        Copyright (C) 2008    Lawrence Berkeley Laboratory
// Author List:
//      Dave Brown 27 May 2011
//------------------------------------------------------------------------
#include "KalmanTests/inc/DetStrawGasElem.hh"
#include "KalmanTests/inc/DetStrawHitType.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "BTrk/TrkBase/TrkDifTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/DetectorModel/DetIntersection.hh"
#include "BTrk/MatEnv/MatDBInfo.hh"
#include <assert.h>
#include <iostream>

namespace mu2e {
  DetStrawGasElem::DetStrawGasElem(DetStrawHitType* stype, TrkStrawHit* strawhit, std::string name) :
    DetElem(stype,name.c_str(),strawhit->straw().index().asInt()),_strawhit(strawhit) {}
  
  DetStrawGasElem::~DetStrawGasElem(){}

  bool
  DetStrawGasElem::reIntersect(const Trajectory* traj,DetIntersection& dinter) const {
    bool retval(false);
    double pdist(0.0);
    double pflt(0.0);
// check if the associated hot is active; if not, use POCA to determine if the track went through the straw
    if(_strawhit->isActive()){
// update the POCA on the hot if necessary.  Must cast the trajectory, FIXME!!!!
      const TrkDifTraj* dtraj = dynamic_cast<const TrkDifTraj*>(traj);
      assert(dtraj != 0);
      _strawhit->updatePoca(dtraj);
      if(_strawhit->pocaStatus().success()){
	// offset the position WRT the hit a little to avoid overlap in the Kalman fit 
	//      double newpath = _strawhit->fltLen() - 0.01*_strawhit->straw().getRadius();
	// check for wild changes in pathlen
	//      static double wildpath = 1000.0;
	//      if(dinter.delem != 0 && fabs(newpath -dinter.pathlen) > wildpath) {
	//	std::cout << "wild path length change in reintersect, old = " << dinter.pathlen << " new = " << newpath << std::endl;
	//	return false;
	//      }
	pflt = _strawhit->fltLen();
	pdist = fabs(_strawhit->driftRadius());
	retval = true;
      }
    } else {
// use POCA to determine if the track went through the straw, even if it's inactive
      TrkPoca poca(*traj,_strawhit->fltLen(),*_strawhit->hitTraj(),_strawhit->hitLen());
      if(poca.status().success() && fabs(poca.doca()) < _strawhit->straw().getRadius()){
	pflt = poca.flt1();
	pdist = fabs(poca.doca());
 	retval = true;
      }
    }
    if(retval){
      // call the base class function
      DetElem::reIntersect(traj,dinter);
      dinter.pathlen = pflt - 0.01*_strawhit->straw().getRadius();
      CLHEP::Hep3Vector tdir = traj->direction(dinter.pathlen);
      double dpath = _strawhit->gasPath(pdist,tdir);

      dinter.pathrange[0] = dinter.pathlen-dpath;
      dinter.pathrange[1] = dinter.pathlen+dpath;
    }
    return retval;
  }
}
