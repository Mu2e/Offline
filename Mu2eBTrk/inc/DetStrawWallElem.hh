//--------------------------------------------------------------------------
// Name:
//   DetStrawWallElem: Dummy class to represent a straw element.
//        Copyright (C) 2008    Lawrence Berkeley Laboratory
// Author List:
//      Dave Brown 27 May 2011
//------------------------------------------------------------------------
#ifndef DetStrawWallElem_hh
#define DetStrawWallElem_hh

#include "BTrk/DetectorModel/DetElem.hh"
#include <string>

namespace mu2e {
  class TrkStrawHit;
  class DetStrawHitType;
// element class  
  class DetStrawWallElem : public DetElem {
  public:
// construct from a TrkStrawHit
    DetStrawWallElem(DetStrawHitType* stype, TrkStrawHit* strawhit, std::string name="DetStrawWallElem");
    virtual ~DetStrawWallElem();
// DetElem interface
    virtual int intersect(const Trajectory*,DetIntersection&) const{ return 0; }
    virtual bool reIntersect(const Trajectory* traj,DetIntersection& dinter) const;
  protected:
    virtual HepPoint coordToPoint( const TypeCoord* aCoord ) const { return HepPoint(0.0,0.0,0.0); }
//  private:
    mutable TrkStrawHit* _strawhit;
  };
}

#endif
