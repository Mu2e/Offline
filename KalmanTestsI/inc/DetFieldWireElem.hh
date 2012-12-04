//
// $Id: DetFieldWireElem.hh,v 1.2 2012/12/04 00:51:26 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:26 $
//
// Original author G. Tassielli
//

#ifndef DetFieldWireElem_hh
#define DetFieldWireElem_hh

#include "KalmanTests/inc/DetStrawWallElem.hh"

namespace mu2e {
  class TrkCellHit;
  class DetStrawHitType;
// element class  
  class DetFieldWireElem : public DetStrawWallElem {
  public:

    enum wireSide { bottom=0, side, top};

// construct from a TrkCellHit
    DetFieldWireElem(DetStrawHitType* stype, TrkCellHit* cellhit, int sideFlag=0, std::string name="DetFieldWireElem");
    virtual ~DetFieldWireElem();
// DetElem interface
//    virtual int intersect(const Trajectory*,DetIntersection&) const{ return 0; }
    virtual bool reIntersect(const Trajectory* traj,DetIntersection& dinter) const;
//  protected:
//    virtual HepPoint coordToPoint( const TypeCoord* aCoord ) const { return HepPoint(0.0,0.0,0.0); }
//  private:
//    mutable TrkStrawHit* _strawhit;

    wireSide _wSide;

  };
}

#endif
