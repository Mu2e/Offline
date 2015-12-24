//--------------------------------------------------------------------------
// Name:
//   DetStrawGasElem: class to represent the straw gas material
//        Copyright (C) 2008    Lawrence Berkeley Laboratory
// Author List:
//      Dave Brown 27 May 2011
//------------------------------------------------------------------------
#ifndef DetStrawGasElem_hh
#define DetStrawGasElem_hh

#include "BTrk/DetectorModel/DetElem.hh"
#include "Mu2eBTrk/inc/DetUniformType.hh"
#include <string>

namespace mu2e {
  class Straw;
  class DetUniformType;
// element class  
  class DetStrawGasElem : public DetElem {
  public:
// construct from a Straw
    DetStrawGasElem(DetUniformType*, Straw const* straw, std::string name);
    virtual ~DetStrawGasElem();
// accessors
    const Straw* straw() const { return _straw; }
// DetElem interface
    virtual int intersect(const Trajectory*,DetIntersection&) const{ return 0; }
    virtual bool reIntersect(const Trajectory* traj,DetIntersection& dinter) const;
    double gasPath(double pdist,CLHEP::Hep3Vector const& tdir) const; // track pathlength through 1/2 the gas of the straw
  protected:
    virtual HepPoint coordToPoint( const TypeCoord* aCoord ) const { return HepPoint(0.0,0.0,0.0); }
//  private:
    const Straw* _straw;
    Trajectory* _wtraj; // wire trajectory
    double _pathoffset; // offset to path to avoid confusion with hits, etc
  };
}

#endif
