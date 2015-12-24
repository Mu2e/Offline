//--------------------------------------------------------------------------
// Name:
//   DetStrawGasElem: class to represent the straw gas material
//        Copyright (C) 2008    Lawrence Berkeley Laboratory
// Author List:
//      Dave Brown 27 May 2011
//------------------------------------------------------------------------
#include "Mu2eBTrk/inc/DetStrawGasElem.hh"
#include "Mu2eBTrk/inc/DetUniformType.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/DetectorModel/DetIntersection.hh"
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include <assert.h>
#include <iostream>

using CLHEP::Hep3Vector;

namespace mu2e {
  DetStrawGasElem::DetStrawGasElem(DetUniformType* stype,const Straw* straw, std::string name) :
    DetElem(stype,name.c_str(),straw->index().asInt()),_straw(straw), _wtraj(0),
    _pathoffset(-0.1) {
// rebuilding the wire traj for each element and the hit is redundant, the traj should be a propoerty of the straw, FIXME!!!
    Hep3Vector const& wiredir = _straw->getDirection();
    Hep3Vector const& mid = _straw->getMidPoint();
    _wtraj = new TrkLineTraj(HepPoint(mid.x(),mid.y(),mid.z()),wiredir,0);
  }

  DetStrawGasElem::~DetStrawGasElem() {
    delete _wtraj;
  }

  bool
  DetStrawGasElem::reIntersect(const Trajectory* traj,DetIntersection& dinter) const {
    bool retval(false);
    double pdist(0.0);
    double pflt(0.0);
// use POCA to determine if the track went through the straw
    TrkPoca poca(*traj,dinter.pathlen,*_wtraj,0);
    if(poca.status().success() && fabs(poca.doca()) < _straw->getRadius()){
      pflt = poca.flt1();
      pdist = fabs(poca.doca());
      retval = true;
      // call the base class function
      DetElem::reIntersect(traj,dinter);
      dinter.pathlen = pflt - 0.01*_straw->getRadius();
      CLHEP::Hep3Vector tdir = traj->direction(dinter.pathlen);
      double dpath = gasPath(pdist,tdir);
      dinter.pathrange[0] = dinter.pathlen-dpath;
      dinter.pathrange[1] = dinter.pathlen+dpath;
    }
    return retval;
  }
  // compute the pathlength through half the gas , given the drift distance and straw geometry
  double DetStrawGasElem::gasPath(double pdist,Hep3Vector const& tdir) const {
    double thick = _straw->getThickness();
    double radius = _straw->getRadius();
    radius-=thick;
    double hlen = _straw->getHalfLength();
    if(pdist>=radius)
      pdist = 0.96*radius;
    double gaspath = sqrt( (radius+pdist)*(radius-pdist) );
// scale for the other dimension
    double cost = tdir.dot(_straw->getDirection());
    if(fabs(cost)<0.999)
      gaspath /= sqrt( (1.0-cost)*(1.0+cost) );
    else
      gaspath = hlen;
// use half-length as maximum length
    gaspath = std::min(gaspath,hlen);
//NAN test
    if(gaspath != gaspath){
      std::cout << "NAN gas" << std::endl;
      gaspath = 2*radius;
    }
    return gaspath;
  }

}

