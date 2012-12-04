//
// $Id: DetSenseWireElem.hh,v 1.2 2012/12/04 00:51:26 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:26 $
//
// Original author G. Tassielli
//

#ifndef DetSenseWireElem_hh
#define DetSenseWireElem_hh

#include "KalmanTests/inc/DetStrawWallElem.hh"

namespace mu2e {
  class TrkCellHit;
  class DetStrawHitType;
// element class  
  class DetSenseWireElem : public DetStrawWallElem {
  public:
// construct from a TrkCellHit
    DetSenseWireElem(DetStrawHitType* stype, TrkCellHit* cellhit, std::string name="DetSenseWireElem");
    virtual ~DetSenseWireElem();
// DetElem interface
//    virtual int intersect(const Trajectory*,DetIntersection&) const{ return 0; }
    virtual bool reIntersect(const Trajectory* traj,DetIntersection& dinter) const;
//  protected:
//    virtual HepPoint coordToPoint( const TypeCoord* aCoord ) const { return HepPoint(0.0,0.0,0.0); }
//  private:
//    mutable TrkStrawHit* _strawhit;

  };
}

#endif
