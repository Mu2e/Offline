//--------------------------------------------------------------------------
// Name:
//   DetStrawWallElem: Dummy class to represent a straw element.
//        Copyright (C) 2008    Lawrence Berkeley Laboratory
// Author List:
//      Dave Brown 27 May 2011
//------------------------------------------------------------------------
#include "KalmanTests/inc/DetStrawWallElem.hh"
#include "KalmanTests/inc/DetStrawHitType.hh"
#include "KalmanTests/inc/TrkStrawHit.hh"
#include "MatEnv/MatDBInfo.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "DetectorModel/DetIntersection.hh"
#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/TrkDifTraj.hh"
#include <assert.h>

namespace mu2e {
//
  DetStrawWallElem::DetStrawWallElem(DetStrawHitType* stype,TrkStrawHit* strawhit) :
  DetElem(stype, "DetStrawWallElem",strawhit->straw().index().asInt()),_strawhit(strawhit) {}
  
  DetStrawWallElem::~DetStrawWallElem(){}

  bool
  DetStrawWallElem::reIntersect(const Trajectory* traj,DetIntersection& dinter) const {
// check if the associated hot is active; if not, return false
    if(!_strawhit->isActive())return false;
// update the POCA on the hot if necessary; must cast the trajectory, FIXME!!!!
    const TrkDifTraj* dtraj = dynamic_cast<const TrkDifTraj*>(traj);
    assert(dtraj != 0);
    _strawhit->updatePoca(dtraj);
// call the base class function
    DetElem::reIntersect(traj,dinter);
// Use the straw hit to define how much material has been traversed.  As the pathlen is
// used to place the material on a Kalman fit, make sure it's not at the identical place as the
// hit, as that causes confusion
    dinter.pathlen = _strawhit->fltLen() + 0.01*_strawhit->straw().getRadius();
// find the trajectory direction at this pathlen
    Hep3Vector tdir = dtraj->direction(_strawhit->fltLen());
    double dpath = _strawhit->wallPath(tdir);
    dinter.pathrange[0] = dinter.pathlen-dpath;
    dinter.pathrange[1] = dinter.pathlen+dpath;
    return true;
  }
}
