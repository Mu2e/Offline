//
// $Id: DetCellGasElem.cc,v 1.2 2012/12/04 00:51:27 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:27 $
//
// Original author G. Tassielli
//

#include "KalmanTestsI/inc/DetCellGasElem.hh"
#include "KalmanTestsI/inc/TrkCellHit.hh"

#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/HelixTraj.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "TrkBase/TrkPoca.hh"
#include "DetectorModel/DetIntersection.hh"
#include "MatEnv/MatDBInfo.hh"

//#include <iostream>

namespace mu2e {
  DetCellGasElem::DetCellGasElem(DetStrawHitType* stype, TrkCellHit* cellhit, std::string name) :
        DetStrawGasElem(stype, (TrkStrawHit*)cellhit, name.c_str()) {}
  
  DetCellGasElem::~DetCellGasElem(){}

  bool
  DetCellGasElem::reIntersect(const Trajectory* traj,DetIntersection& dinter) const {
    bool retval(false);
    //double pdist(0.0);
    double pflt(0.0);
    TrkCellHit *cellhit = (TrkCellHit *)_strawhit;
// check if the associated hot is active; if not, use POCA to determine if the track went through the cell
    if(cellhit->isActive()){
// update the POCA on the hot if necessary.  Must cast the trajectory, FIXME!!!!
      const TrkDifTraj* dtraj = dynamic_cast<const TrkDifTraj*>(traj);
      assert(dtraj != 0);
      cellhit->updatePoca(dtraj);
      if(cellhit->pocaStatus().success()){
	// offset the position WRT the hit a little to avoid overlap in the Kalman fit 
	//      double newpath = cellhit->fltLen() - 0.01*cellhit->straw().getRadius();
	// check for wild changes in pathlen
	//      static double wildpath = 1000.0;
	//      if(dinter.delem != 0 && fabs(newpath -dinter.pathlen) > wildpath) {
	//	std::cout << "wild path length change in reintersect, old = " << dinter.pathlen << " new = " << newpath << std::endl;
	//	return false;
	//      }
	pflt = cellhit->fltLen();
	//pdist = fabs(cellhit->driftRadius());
	retval = true;
      }
    } else {
// use POCA to determine if the track went through the cell, even if it's inactive
      TrkPoca poca(*traj,cellhit->fltLen(),*cellhit->hitTraj(),cellhit->hitLen());
      if(poca.status().success() && fabs(poca.doca()) < cellhit->straw().getRadius()){
	pflt = poca.flt1();
	//pdist = fabs(poca.doca());
 	retval = true;
      }
    }
    if(retval){
      // call the base class function
      DetElem::reIntersect(traj,dinter);
      dinter.pathlen = pflt;// - 0.01*cellhit->straw().getRadius();
      HelixTraj* htraj = ((HelixTraj*) traj);
      Hep3Vector tdir = htraj->direction(dinter.pathlen);
      const HepPoint &tpos = htraj->position(dinter.pathlen);

      /*double dpath =*/ cellhit->cellGasPath(dinter.pathlen,tdir,tpos/*,htraj->omega(),htraj->cosDip(),htraj->angle(dinter.pathlen)*/);
      dinter.pathrange[0] = dinter.pathlen+cellhit->entryDeltaPath();
      dinter.pathrange[1] = dinter.pathlen+cellhit->exitDeltaPath();
      //std::cout<<"gas path "<< dinter.pathlen <<" min "<< dinter.pathrange[0] <<" max "<<dinter.pathrange[1]<<std::endl;

    }
    return retval;
  }
}
