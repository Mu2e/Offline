//--------------------------------------------------------------------------
// Name:
//   DetStrawWallElem: class to represent a straw wall material
//        Copyright (C) 2008    Lawrence Berkeley Laboratory
// Author List:
//      Dave Brown 27 May 2011
//------------------------------------------------------------------------
#ifndef DetStrawWallElem_hh
#define DetStrawWallElem_hh

#include "BTrk/DetectorModel/DetElem.hh"
#include "Mu2eBTrk/inc/DetUniformType.hh"
#include <string>
class Trajectory;
namespace mu2e {
  class Straw;
  class DetUniformType;
// element class  
  class DetStrawWallElem : public DetElem {
  public:
// construct from a Straw
    DetStrawWallElem(DetUniformType* type, const Straw* straw, std::string name);
    virtual ~DetStrawWallElem();
// accessors
    const Straw* straw() const { return _straw; }
// DetElem interface; only reIntersect used.
    virtual int intersect(const Trajectory*,DetIntersection&) const{ return 0; }
    virtual bool reIntersect(const Trajectory* traj,DetIntersection& dinter) const;
    double wallPath(double pdist,CLHEP::Hep3Vector const& tdir) const; // track pathlength through one wall of the straw
  protected:
    virtual HepPoint coordToPoint( const TypeCoord* aCoord ) const { return HepPoint(0.0,0.0,0.0); }
  private:
    const Straw* _straw;
    Trajectory* _wtraj; // wire trajectory
    double _rtol; // tolerance on being inside the straw
    double _pathoffset; // offset to path to avoid confusion with hits, etc
  };
}

#endif
