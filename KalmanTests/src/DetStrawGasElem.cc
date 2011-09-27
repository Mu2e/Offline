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
#include "TrkBase/TrkDifTraj.hh"
#include "TrkBase/TrkErrCode.hh"
#include "DetectorModel/DetIntersection.hh"
#include "MatEnv/MatDBInfo.hh"
#include <assert.h>

namespace mu2e {
  DetStrawGasElem::DetStrawGasElem(DetStrawHitType* stype, TrkStrawHit* strawhit) :
    DetElem(stype,"DetStrawGasElem",strawhit->straw().index().asInt()),_strawhit(strawhit) {}
  
  DetStrawGasElem::~DetStrawGasElem(){}

  bool
  DetStrawGasElem::reIntersect(const Trajectory* traj,DetIntersection& dinter) const {
// check if the associated hot is active; if not, return false
    if(!_strawhit->isActive())return false;
// update the POCA on the hot if necessary.  Must cast the trajectory, FIXME!!!!
    const TrkDifTraj* dtraj = dynamic_cast<const TrkDifTraj*>(traj);
    assert(dtraj != 0);
    _strawhit->updatePoca(dtraj);
// call the base class function
    DetElem::reIntersect(traj,dinter);
// offset the position WRT the hit a little to avoid overlap in the Kalman fit 
    dinter.pathlen = _strawhit->fltLen() - 0.01*_strawhit->straw().getRadius();
// Use the straw hit to define how much material has been traversed
    Hep3Vector tdir = dtraj->direction(_strawhit->fltLen());
    double dpath = _strawhit->gasPath(tdir);
    dinter.pathrange[0] = dinter.pathlen-dpath;
    dinter.pathrange[1] = dinter.pathlen+dpath;
    return true;
  }
}
