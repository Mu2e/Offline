//
// $Id: DetFieldWireElem.cc,v 1.2 2012/12/04 00:51:27 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:27 $
//
// Original author G. Tassielli
//

#include "KalmanTestsI/inc/DetFieldWireElem.hh"
#include "KalmanTestsI/inc/TrkCellHit.hh"

#include "TrkBase/TrkErrCode.hh"
#include "TrkBase/HelixTraj.hh"
#include "TrkBase/TrkDifTraj.hh"
#include "TrkBase/TrkPoca.hh"
#include "DetectorModel/DetIntersection.hh"
#include "MatEnv/MatDBInfo.hh"

//#include <iostream>
#include "cetlib/pow.h"

namespace mu2e {
  DetFieldWireElem::DetFieldWireElem(DetStrawHitType* stype,TrkCellHit* cellhit, int sideFlag, std::string name) :
                DetStrawWallElem(stype, (TrkStrawHit*)cellhit, name.c_str()) {
          if (sideFlag==1) {
                  _wSide = side;
          } else if (sideFlag==2) {
                  _wSide = top;
          } else {
                  _wSide = bottom;
          }
  }
  
  DetFieldWireElem::~DetFieldWireElem(){}

  bool
  DetFieldWireElem::reIntersect(const Trajectory* traj,DetIntersection& dinter) const {
    bool retval(false);
    //double pdist(0.0);
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
	//pdist = fabs(cellhit->poca()->doca());
	retval =  true;
      }
    } else {
// use POCA to determine if the track went through the cell, even if it's inactive
      TrkPoca poca(*traj,cellhit->fltLen(),*cellhit->hitTraj(),cellhit->hitLen());
      if(poca.status().success() && fabs(poca.doca()) < cellhit->cellHandle()->GetCellRad()/*cellhit->straw().getRadius()*/){
	pflt = poca.flt1();
	//pdist = fabs(poca.doca());
	retval = true;
      }
    }
    if(retval){
      // call the base class function
      // As the pathlen is
      // used to place the material on a Kalman fit, make sure it's not at the identical place as the
      // hit, as that causes confusion
      DetElem::reIntersect(traj,dinter);
      dinter.pathlen = pflt;
      HelixTraj* htraj = ((HelixTraj*) traj);
      Hep3Vector tdir = htraj->direction(dinter.pathlen);
      const HepPoint &tpos = htraj->position(dinter.pathlen);
      double dpath = cellhit->fieldWirePath(dinter.pathlen,tdir,tpos/*,htraj->omega(),htraj->cosDip(),htraj->angle(dinter.pathlen)*/);

      if (!cellhit->pcaOutOfCut()) {
              const HepPoint &outEntPos = htraj->position(dinter.pathlen+cellhit->entryDeltaPath());
              HepGeom::Point3D<double> tmpGlobalPos(outEntPos.x(),outEntPos.y(),outEntPos.z());
              HepGeom::Point3D<double> tmpLocalPos;
              cellhit->cellHandle()->Global2Local(tmpGlobalPos,tmpLocalPos);
              const HepPoint &outExtPos = htraj->position(dinter.pathlen+cellhit->exitDeltaPath());
              tmpGlobalPos.set(outExtPos.x(),outExtPos.y(),outExtPos.z());
              cellhit->cellHandle()->Global2Local(tmpGlobalPos,tmpLocalPos);
              cellhit->cellHandle()->Local2Global(tmpLocalPos,tmpGlobalPos);

              double tolerance=0.002;

              CLHEP::Hep3Vector fwPca;
              if ( _wSide==bottom ) {
                      if (cellhit->isGrowingRad()) {dinter.pathlen+=cellhit->entryDeltaPath();}
                      else  {dinter.pathlen+=cellhit->exitDeltaPath();}
                      const HepPoint &tmpEntPos = htraj->position(dinter.pathlen);
                      const CLHEP::Hep3Vector entPos(tmpEntPos.x(),tmpEntPos.y(),tmpEntPos.z());
                      /*CLHEP::Hep3Vector entPos;
                      cellhit->cellHandle()->WirePosAtZ(tmpEntPos.z(),entPos);
                      std::cout<<"entPos on center wire "<<entPos<<std::endl;
                      double rsw = sqrt(cet::sum_of_squares(tmpEntPos.x(),tmpEntPos.y()));
                      entPos.setX( (rsw-cellhit->cellHandle()->GetCellInsideRad())*entPos.x()/rsw );
                      entPos.setY( (rsw-cellhit->cellHandle()->GetCellInsideRad())*entPos.y()/rsw );
                      std::cout<<"entPos moved "<<entPos<<std::endl;*/

                      Hep3Vector entDir = htraj->direction(dinter.pathlen);
                      CellGeometryHandle::FWireSide sideFlag(CellGeometryHandle::bottom);
                      dpath = cellhit->cellHandle()->CrossingPathOnFieldWires(entPos,entDir,sideFlag,fwPca,tolerance);
                      cellhit->_hitFWire_b=dpath>0.0;
                      cellhit->_pathlInFWire_b=dpath;
                      if(cellhit->_hitFWire_b) {
                              std::cout<<"Fw wire hit for cell "<<cellhit->cellHandle()->GetITCell()->Id()<<std::endl;
                              std::cout<<"side fwire "<<_wSide<<" pathlen "<<dinter.pathlen<<" dpath "<<dpath<<std::endl;
                      }
              } else if ( _wSide==side ) {
                      dinter.pathlen+=cellhit->entryDeltaPath();
                      const HepPoint &tmpEntPos = htraj->position(dinter.pathlen);
                      const CLHEP::Hep3Vector entPos(tmpEntPos.x(),tmpEntPos.y(),tmpEntPos.z());
                      /*CLHEP::Hep3Vector entPos;
                      cellhit->cellHandle()->WirePosAtZ(tmpEntPos.z(),entPos);
                      entPos.setX(entPos.x()-cellhit->cellHandle()->GetCellInsideRad());*/

                      Hep3Vector entDir = htraj->direction(dinter.pathlen);
                      CellGeometryHandle::FWireSide sideFlag(CellGeometryHandle::side);
                      dpath = cellhit->cellHandle()->CrossingPathOnFieldWires(entPos,entDir,sideFlag,fwPca,tolerance);
                      cellhit->_hitFWire_s=dpath>0.0;
                      cellhit->_pathlInFWire_s=dpath;
                      if(cellhit->_hitFWire_s) {
                              std::cout<<"Fw wire hit for cell "<<cellhit->cellHandle()->GetITCell()->Id()<<std::endl;
                              std::cout<<"side fwire "<<_wSide<<" pathlen "<<dinter.pathlen<<" dpath "<<dpath<<std::endl;
                      }
              } else if ( _wSide==top ){
                      if (cellhit->isGrowingRad()) {dinter.pathlen+=cellhit->exitDeltaPath();}
                      else  {dinter.pathlen+=cellhit->entryDeltaPath();}
                      const HepPoint &tmpEntPos = htraj->position(dinter.pathlen);
                      const CLHEP::Hep3Vector entPos(tmpEntPos.x(),tmpEntPos.y(),tmpEntPos.z());
                      /*CLHEP::Hep3Vector entPos;
                      cellhit->cellHandle()->WirePosAtZ(tmpEntPos.z(),entPos);
                      double rsw = sqrt(cet::sum_of_squares(tmpEntPos.x(),tmpEntPos.y()));
                      entPos.setX( (rsw+cellhit->cellHandle()->GetCellInsideRad())*entPos.x()/rsw );
                      entPos.setY( (rsw+cellhit->cellHandle()->GetCellInsideRad())*entPos.y()/rsw );*/

                      Hep3Vector entDir = htraj->direction(dinter.pathlen);
                      CellGeometryHandle::FWireSide sideFlag(CellGeometryHandle::top);
                      dpath = cellhit->cellHandle()->CrossingPathOnFieldWires(entPos,entDir,sideFlag,fwPca,tolerance);
                      cellhit->_hitFWire_t=dpath>0.0;
                      cellhit->_pathlInFWire_t=dpath;
                      if(cellhit->_hitFWire_t) {
                              std::cout<<"Fw wire hit for cell "<<cellhit->cellHandle()->GetITCell()->Id()<<std::endl;
                              std::cout<<"side fwire "<<_wSide<<" pathlen "<<dinter.pathlen<<" dpath "<<dpath<<std::endl;
                      }
              }
              dpath*=0.5;
              dinter.pathrange[0] = dinter.pathlen-dpath;
              dinter.pathrange[1] = dinter.pathlen+dpath;
      } else {
              dinter.pathrange[0] = dinter.pathlen;
              dinter.pathrange[1] = dinter.pathlen;
      }


    }
    return retval;
  }
}

