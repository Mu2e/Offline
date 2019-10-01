//--------------------------------------------------------------------------
// Name:
//   DetStrawElem: class to represent the straw gas material
//        Copyright (C) 2008    Lawrence Berkeley Laboratory
// Author List:
//      Dave Brown 27 May 2011
//------------------------------------------------------------------------
#include "Mu2eBTrk/inc/DetStrawElem.hh"
#include "Mu2eBTrk/inc/DetStrawType.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "BTrk/TrkBase/TrkHit.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/DetectorModel/DetIntersection.hh"
#include "BTrk/BbrGeom/TrkLineTraj.hh"
#include "cetlib_except/coded_exception.h"
#include <assert.h>
#include <iostream>

using CLHEP::Hep3Vector;

namespace mu2e {
  DetStrawElem::DetStrawElem(DetStrawType* stype,const Straw* straw) :
    DetElem(stype,"TrackerStraw",straw->id().asUint16()),_straw(straw), _wtraj(0), _stype(stype){
// the traj should be a propoerty of the straw, FIXME!!!
    Hep3Vector const& wiredir = _straw->getDirection();
    Hep3Vector const& mid = _straw->getMidPoint();
    _wtraj = new TrkLineTraj(HepPoint(mid.x(),mid.y(),mid.z()),wiredir,0);
  }

  DetStrawElem::~DetStrawElem() {
    delete _wtraj;
  }

  int
  DetStrawElem::intersect(const Trajectory* traj, DetIntersection& dinter) const {
    if(reIntersect(traj,dinter))
      return 1;
    else
      return 0;
  }

  bool
  DetStrawElem::reIntersect(const Trajectory* traj,DetIntersection& dinter) const {
    bool retval(false);
// use POCA to determine if the track went through the straw.  The DetIntersection
// tolerance defines how close we have to get
    TrkPoca poca(*traj,dinter.pathlen,*_wtraj,0);
    if(poca.status().success()){
// if there's an associated active hit with this element, then force the intersection
  // otherwise, require the distance to the wire be inside the straw (within tolerance)
  // should also check the distance along the wire is inside the active length FIXME!!
      if( (dinter.thit != 0 && dinter.thit->isActive()) ||
	fabs(poca.doca()) < _straw->getRadius() + _stype->tolerance() ){
	retval = true;
// call the base class function to update the base members
	DetElem::reIntersect(traj,dinter);
// fill the DetIntersection.  Use the gas path give the intersection for the whole straw
	CLHEP::Hep3Vector tdir = traj->direction(dinter.pathlen);
	double dpath = gasPath(fabs(poca.doca()),tdir);
	dinter.pathlen = poca.flt1() + _stype->offset();
	dinter.dist = poca.doca(); // NB: this can be negative!
	dinter.pathrange[0] = dinter.pathlen-dpath;
	dinter.pathrange[1] = dinter.pathlen+dpath;
      }
    }
    return retval;
  }

// compute the material effects of traversing the straw.  This includes
// the straw gas and walls (and eventually wire, FIXME!!)
  void DetStrawElem::materialInfo(const DetIntersection& dinter,
      double momentum,
      TrkParticle const& tpart,
      double& deflectRMS,
      double& pfracRMS,
      double& pfrac,
      trkDirection dedxdir) const {
// compute the path through the straw wall and gas (and eventually test for wire intersections!)
    CLHEP::Hep3Vector tdir = dinter.trajet->direction(dinter.pathlen);
    double gaspath = gasPath(dinter.dist,tdir);
    double wallpath = wallPath(dinter.dist,tdir);
// compute the material info for these materials using the base class function
    double gasdeflectRMS, gaspfracRMS,gaspfrac;
    DetElem::materialInfo(*_stype->gasMaterial(),2*gaspath,momentum,tpart,gasdeflectRMS,gaspfracRMS,gaspfrac,dedxdir);
    double walldeflectRMS, wallpfracRMS,wallpfrac;
    DetElem::materialInfo(*_stype->wallMaterial(),2*wallpath,momentum,tpart,walldeflectRMS,wallpfracRMS,wallpfrac,dedxdir);
  // combine these to give the aggregate effect
    deflectRMS = sqrt(gasdeflectRMS*gasdeflectRMS + walldeflectRMS*walldeflectRMS);
    pfracRMS = sqrt(gaspfracRMS*gaspfracRMS + wallpfracRMS*wallpfracRMS);
    pfrac = gaspfrac + wallpfrac;
  }

  // compute the pathlength through half the gas , given the drift distance and straw geometry
  double DetStrawElem::gasPath(double pdist,Hep3Vector const& tdir) const {
    double radius = _straw->getRadius();
    double hlen = _straw->halfLength();
// if the POCA distance is outside or too close the outside of the straw, force it inside
    pdist = std::min(fabs(pdist),_stype->maxRadiusFraction()*radius);
    double gaspath = sqrt( (radius+pdist)*(radius-pdist) );
// scale for the other dimension.  Maximum path is a fraction of the straw length
    double cost = tdir.dot(_straw->getDirection());
// avoid degenerate case
    if(fabs(cost)<0.999)
      gaspath /= sqrt( (1.0-cost)*(1.0+cost) );
    else
      gaspath = hlen;
// restrict to sensible physical distance
    gaspath = std::min(gaspath,hlen);
    return gaspath;
  }

  double
  DetStrawElem::radiationFraction(const DetIntersection& dinter) const {
// compute the path through the straw wall and gas (and eventually test for wire intersections!)
    CLHEP::Hep3Vector tdir = dinter.trajet->direction(dinter.pathlen);
    double gaspath = gasPath(dinter.dist,tdir);
    double wallpath = wallPath(dinter.dist,tdir);
    double retval = _stype->gasMaterial()->radiationFraction(2*gaspath);
    retval += _stype->wallMaterial()->radiationFraction(2*wallpath);
    return retval;
  }

// compute the pathlength through one wall of the straw, given the drift distance and straw geometry
  double DetStrawElem::wallPath(double pdist,Hep3Vector const& tdir) const {
    double thick = _straw->getThickness();
    double radius = _straw->getRadius();
    double inRadius = radius - thick;
// if the POCA distance is outside or too close the outside of the straw, force it inside
    pdist = std::min(fabs(pdist),_stype->maxRadiusFraction()*inRadius);
    double wallpath =  (sqrt( (radius+pdist)*(radius-pdist) ) -
	sqrt( (inRadius+pdist)*(inRadius-pdist) ));
    // scale for the other dimension
    double cost = tdir.dot(_straw->getDirection());
// avoid degenerate case
    if(fabs(cost)<0.999)
      wallpath /= sqrt( (1.0-cost)*(1.0+cost) );
    else
      wallpath = radius;
  // restrict to a sensible maximu distance
    wallpath = std::min(wallpath,radius);
    return wallpath;
  }
}
