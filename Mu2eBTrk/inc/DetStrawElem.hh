//--------------------------------------------------------------------------
// Name:
//   DetStrawElem: class to represent a straw. Depending on the 
//   type this class can represent the material of the gas, wall, wire, or all
//        Copyright (C) 2008    Lawrence Berkeley Laboratory
// Author List:
//      Dave Brown 27 May 2011
//------------------------------------------------------------------------
#ifndef DetStrawElem_hh
#define DetStrawElem_hh

#include "BTrk/DetectorModel/DetElem.hh"
#include <string>

namespace mu2e {
  class Straw;
  class DetStrawType;
// element class  
  class DetStrawElem : public DetElem {
  public:
// construct from a Straw
    DetStrawElem(DetStrawType* strawtype, Straw const* straw);
    virtual ~DetStrawElem();
// DetElem interface; intersect just calls down to reintersect.  Note that the entrance
// flag of the DetIntersection defines what parts of the material to model
    virtual int intersect(const Trajectory*,DetIntersection&) const;
    virtual bool reIntersect(const Trajectory* traj,DetIntersection& dinter) const;
// overwrite the material info function to allow heterogenous materials with a single intersection
    virtual void materialInfo(const DetIntersection&,
	double momentum,
	TrkParticle const& tpart,
	double& deflectRMS,
	double& pFractionRMS,
	double& pFraction,
	trkDirection dedxdir=trkOut) const;
    // override KalMaterial function which is broken and shouldn't be used
    double radiationFraction(const DetIntersection&) const;
    // accessors
    const Straw* straw() const { return _straw; }
    const Trajectory* wireTraj() const { return _wtraj; }
// specific functions for different material intersections
    double gasPath(double pdist,CLHEP::Hep3Vector const& tdir) const; // track pathlength through 1/2 the gas of the straw, given distance to the wire center
    double wallPath(double pdist,CLHEP::Hep3Vector const& tdir) const; // track pathlength through one wall of the straw, given distance to the wire center
  private:
    virtual HepPoint coordToPoint( const TypeCoord* aCoord ) const { return HepPoint(0.0,0.0,0.0); }
    const Straw* _straw;
    Trajectory* _wtraj; // wire trajectory, owned by this object;  This should be in the straw, FIXME!
    // remember the fully-typed type
    const DetStrawType* _stype;
  };
}

#endif
