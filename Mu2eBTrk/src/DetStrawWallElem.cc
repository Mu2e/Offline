//--------------------------------------------------------------------------
// Name:
//   DetStrawWallElem: class to represent a straw wall material
//        Copyright (C) 2008    Lawrence Berkeley Laboratory
// Author List:
//      Dave Brown 27 May 2011
//------------------------------------------------------------------------
#include "Mu2eBTrk/inc/DetStrawWallElem.hh"
#include "Mu2eBTrk/inc/DetUniformType.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/DetectorModel/DetIntersection.hh"
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include <assert.h>

using CLHEP::Hep3Vector;

namespace mu2e {
//
  DetStrawWallElem::DetStrawWallElem(DetUniformType* stype,const Straw* straw, std::string name) :
  DetElem(stype, name.c_str(),straw->index().asInt()),_straw(straw),_wtraj(0),
  _rtol(0.1), _pathoffset(0.1) {
  // rebuilding the wire traj for each element and the hit is redundant, the traj should be a propoerty of the straw, FIXME!!!
    Hep3Vector const& wiredir = _straw->getDirection();
    Hep3Vector const& mid = _straw->getMidPoint();
    _wtraj = new TrkLineTraj(HepPoint(mid.x(),mid.y(),mid.z()),wiredir,0);
  }
  
  DetStrawWallElem::~DetStrawWallElem(){
    delete _wtraj;
  }

  bool
  DetStrawWallElem::reIntersect(const Trajectory* traj,DetIntersection& dinter) const {
    bool retval(false);
// use POCA to see how far from the wire the track went.
    TrkPoca poca(*traj,dinter.pathlen,*_wtraj,0);
    if(poca.status().success() && fabs(poca.doca()) < _straw->getRadius() + _rtol){
      double pflt = poca.flt1();
      double pdist = fabs(poca.doca());
      retval = true;
      // used to place the material on a Kalman fit, make sure it's not at the identical place as the
      // hit, as that causes confusion
      DetElem::reIntersect(traj,dinter);
      dinter.pathlen = pflt + _pathoffset;
      CLHEP::Hep3Vector tdir = traj->direction(dinter.pathlen);
      double dpath = wallPath(pdist,tdir);
      dinter.pathrange[0] = dinter.pathlen-dpath;
      dinter.pathrange[1] = dinter.pathlen+dpath;
    }
    return retval;
  }
// compute the pathlength through one wall of the straw, given the drift distance and straw geometry
  double DetStrawWallElem::wallPath(double pdist,Hep3Vector const& tdir) const {
    double thick = _straw->getThickness();
    double radius = _straw->getRadius();
    double inRadius = radius-thick;
    if(pdist>=inRadius)
      pdist = 0.96*inRadius;
    double wallpath =  (sqrt( (radius+pdist)*(radius-pdist) ) -
	sqrt( (inRadius+pdist)*(inRadius-pdist) ));
    // scale for the other dimension
    double cost = tdir.dot(_straw->getDirection());
    if(fabs(cost)<0.999)
      wallpath /= sqrt( (1.0-cost)*(1.0+cost) );
    else
      wallpath = radius;
    // use half-length as maximum length
    wallpath = std::min(wallpath,radius);
    // test for NAN	
    if(wallpath != wallpath){
      std::cout << "NAN wall" << std::endl;
      wallpath = thick;
    }
    return wallpath;
  }
}

