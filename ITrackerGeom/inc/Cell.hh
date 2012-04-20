#ifndef ITrackerGeom_Cell_hh
#define ITrackerGeom_Cell_hh

#include <deque>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "ITrackerGeom/inc/CellDetail.hh"
#include "ITrackerGeom/inc/CellId.hh"
#include "ITrackerGeom/inc/Wire.hh"

#include "TrackerGeom/inc/Straw.hh"

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

class Cell : public Straw {

  friend class CellGeometryHandle;
  friend class CellGeometryHandle_v2;
  friend class CellGeometryHandle_v3;
  friend class CellGeometryHandle_v2_DBL;
  friend class CellGeometryHandle_v3_DBL;
  friend class ITracker;
  friend class ITrackerMaker;
  friend class SuperLayer;

public:

  // A free function, returning void, that takes a const Cell& as an argument.
  typedef void (*CellFunction)( const Cell& s);

  Cell();

  // Constructor using sense wire info.
  Cell( CellId id,
         boost::shared_ptr<CellDetail> detail,
         boost::shared_ptr<Wire> senseWire
         );

//  // Constructor using sense wire info.
//  Cell( CellId id,
//         CellIndex index,
//         boost::shared_ptr<CellDetail> detail,
//         int detailIndex,
//         boost::shared_ptr<Wire> senseWire
//         );

  ~Cell ();

  CellId Id() const { return _id;}
//  CellIndex Index() const { return _index;}

  boost::shared_ptr<CellDetail> getDetail() const { return _detail;}

  boost::shared_ptr<Wire> getWire() const { return  _senseWire; }

  const CLHEP::Hep3Vector& getMidPoint()  const { return _tmpMidPoint;  /*return _senseWire.get()->getMidPoint();*/}

  const CLHEP::Hep3Vector& getDirection() const { return _tmpDirection; /*return _senseWire.get()->getDirection();*/}

  double getHalfLength() const { return _senseWire.get()->getDetail()->halfLength();}

  double getRadius() const { return getDetail()->CirumscribedRadius();}
  double getThickness() const { return 0.;}
//  int hack;

protected:

  // Identifier
  CellId _id;

  // Detailed description of a cell.
  boost::shared_ptr<CellDetail> _detail;
  int _detailIndex;

  boost::shared_ptr<Wire> _senseWire;

  CLHEP::Hep3Vector _tmpMidPoint;
  CLHEP::Hep3Vector _tmpDirection;

};

}  //namespace mu2e

#endif /* ITrackerGeom_Cell_hh */
