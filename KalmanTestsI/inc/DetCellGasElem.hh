//
// DetCellGasElem: Dummy class to represent a cell element.
//
// $Id: DetCellGasElem.hh,v 1.2 2012/12/04 00:51:26 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:26 $
//
// Original author G. Tassielli
//

#ifndef DetCellGasElem_hh
#define DetCellGasElem_hh

#include "KalmanTests/inc/DetStrawGasElem.hh"

namespace mu2e {
  class TrkCellHit;
  class DetStrawHitType;
// element class  
  class DetCellGasElem : public DetStrawGasElem {
  public:
// construct from a TrkCellHit
    DetCellGasElem(DetStrawHitType*, TrkCellHit* cellhit, std::string name="DetCellGasElem");
    virtual ~DetCellGasElem();
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
