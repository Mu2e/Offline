//
// $Id: DetSenseWireElem.cc,v 1.2 2012/12/04 00:51:27 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:27 $
//
// Original author G. Tassielli
//

#include "KalmanTestsI/inc/DetSenseWireElem.hh"
#include "KalmanTestsI/inc/TrkCellHit.hh"

#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/HelixTraj.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "TrkBase/TrkPoca.hh"
#include "DetectorModel/DetIntersection.hh"
#include "MatEnv/MatDBInfo.hh"

//#include <iostream>

namespace mu2e {
  DetSenseWireElem::DetSenseWireElem(DetStrawHitType* stype,TrkCellHit* cellhit, std::string name) :
                DetStrawWallElem(stype, (TrkStrawHit*)cellhit, name.c_str()) {}
  
  DetSenseWireElem::~DetSenseWireElem(){}

  bool
  DetSenseWireElem::reIntersect(const Trajectory* traj,DetIntersection& dinter) const {
    bool retval(false);
    double pdist(0.0);
    double pflt(0.0);

    TrkCellHit *cellhit = (TrkCellHit *)_strawhit;
// check if the associated hot is active; if not, return false
    if(cellhit->isActive()){
      // update the POCA on the hot if necessary; must cast the trajectory, FIXME!!!!
      const TrkDifTraj* dtraj = dynamic_cast<const TrkDifTraj*>(traj);
      assert(dtraj != 0);
      cellhit->updatePoca(dtraj);
      if(cellhit->pocaStatus().success()){
	pflt = cellhit->fltLen();
// can't use drift radius as t2d not initialized on construction
//	pdist = fabs(cellhit->driftRadius());
	pdist = fabs(cellhit->poca()->doca());
	retval =  true;
      }
    } else {
// use POCA to determine if the track went through the cell, even if it's inactive
      TrkPoca poca(*traj,cellhit->fltLen(),*cellhit->hitTraj(),cellhit->hitLen());
      if(poca.status().success() && fabs(poca.doca()) < cellhit->cellHandle()->GetCellRad()/*cellhit->straw().getRadius()*/){
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
      dinter.pathlen = pflt;// + 0.01*cellhit->cellHandle()->GetCellRad()/*cellhit->straw().getRadius()*/;
      HelixTraj* htraj = ((HelixTraj*) traj);
      Hep3Vector tdir = htraj->direction(dinter.pathlen);
      const HepPoint &tpos = htraj->position(dinter.pathlen);
      double dpath = cellhit->senseWirePath(pdist,tdir,tpos);
      dinter.pathlen+=dpath;
      CLHEP::Hep3Vector swPca;
      const HepPoint &tmpClsPos = htraj->position(dinter.pathlen);
      const CLHEP::Hep3Vector clsPos(tmpClsPos.x(),tmpClsPos.y(),tmpClsPos.z());
      /*CLHEP::Hep3Vector clsPos;
      cellhit->cellHandle()->WirePosAtZ(tmpClsPos.z(),clsPos);*/

      Hep3Vector clsDir = htraj->direction(dinter.pathlen);
      double tolerance=0.002;
      dpath = cellhit->cellHandle()->CrossingPathOnSenseWires(clsPos,clsDir,swPca,tolerance);
      cellhit->_hitSWire=dpath>0.0;
      cellhit->_pathlInSWire=dpath;
      dpath*=0.5;

      dinter.pathrange[0] = dinter.pathlen-dpath;
      dinter.pathrange[1] = dinter.pathlen+dpath;
      if(cellhit->_hitSWire) {
              std::cout<<"Sw wire hit for cell "<<cellhit->cellHandle()->GetITCell()->Id()<<std::endl;
              std::cout<<"sense wire: pathlen "<<dinter.pathlen<<" dpath "<<dpath<<std::endl;
      }
    }
    return retval;
  }
}

