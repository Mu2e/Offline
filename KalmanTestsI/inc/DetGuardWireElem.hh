//
// $Id: DetGuardWireElem.hh,v 1.2 2012/12/04 00:51:26 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/12/04 00:51:26 $
//
// Original author G. Tassielli
//

#ifndef DetGuardWireElem_hh
#define DetGuardWireElem_hh

#include "DetectorModel/DetElem.hh"
#include "ITrackerGeom/inc/Wire.hh"
#include <string>
#include <vector>
#include <utility>

namespace mu2e {
// element class  
  class DetStrawHitType;
  class DetGuardWireElem : public DetElem {
  public:
// construct
    DetGuardWireElem(DetStrawHitType*, int id, double zmin, double zmax, const std::vector<boost::shared_ptr<Wire> > &fieldWires, std::string name="DetGuardWireElem");
    virtual ~DetGuardWireElem();
// DetElem interface
    virtual int intersect(const Trajectory*,DetIntersection&) const{ return 0; }
    virtual bool reIntersect(const Trajectory* traj,DetIntersection& dinter) const;
    //void addWire(CLHEP::Hep3Vector *midP, CLHEP::Hep3Vector *dir) {
    //        _wiresMP_Dir.push_back( std::make_pair<CLHEP::Hep3Vector*>(midP,dir) );
    //}
  protected:
    virtual HepPoint coordToPoint( const TypeCoord* aCoord ) const { return HepPoint(0.0,0.0,0.0); }

    double _zmin;
    double _zmax;
    const std::vector<boost::shared_ptr<Wire> > _wires;

  };
}

#endif
