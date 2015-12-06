//--------------------------------------------------------------------------
// Name:
//   DetStrawWallElem: Dummy class to represent a straw element.
//        Copyright (C) 2008    Lawrence Berkeley Laboratory
// Author List:
//      Dave Brown 27 May 2011
//------------------------------------------------------------------------
#include "Mu2eBTrk/inc/DetStrawWallElem.hh"
#include "Mu2eBTrk/inc/DetStrawHitType.hh"
#include "Mu2eBTrk/inc/TrkStrawHit.hh"
#include "BTrk/MatEnv/MatDBInfo.hh"
#include "BTrk/TrkBase/TrkDifTraj.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/DetectorModel/DetIntersection.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/TrkBase/TrkDifTraj.hh"
#include <assert.h>

namespace mu2e {
//
  DetStrawWallElem::DetStrawWallElem(DetStrawHitType* stype,TrkStrawHit* strawhit, std::string name) :
  DetElem(stype, name.c_str(),strawhit->straw().index().asInt()),_strawhit(strawhit) {}
  
  DetStrawWallElem::~DetStrawWallElem(){}

  bool
  DetStrawWallElem::reIntersect(const Trajectory* traj,DetIntersection& dinter) const {
    bool retval(false);
    double pdist(0.0);
    double pflt(0.0);
// check if the associated hot is active; if not, return false
    if(_strawhit->isActive()){
      // update the POCA on the hot if necessary; must cast the trajectory, FIXME!!!!
      const TrkDifTraj* dtraj = dynamic_cast<const TrkDifTraj*>(traj);
      assert(dtraj != 0);
      _strawhit->updatePoca(dtraj);
      if(_strawhit->poca().status().success()){
	pflt = _strawhit->fltLen();
// can't use drift radius as t2d not initialized on construction
//	pdist = fabs(_strawhit->driftRadius());
	pdist = fabs(_strawhit->poca().doca());
	retval =  true;
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
      // As the pathlen is
      // used to place the material on a Kalman fit, make sure it's not at the identical place as the
      // hit, as that causes confusion
      DetElem::reIntersect(traj,dinter);
      dinter.pathlen = pflt + 0.01*_strawhit->straw().getRadius();
      CLHEP::Hep3Vector tdir = traj->direction(dinter.pathlen);
      double dpath = _strawhit->wallPath(pdist,tdir);
      dinter.pathrange[0] = dinter.pathlen-dpath;
      dinter.pathrange[1] = dinter.pathlen+dpath;
    }
    return retval;
  }
}

