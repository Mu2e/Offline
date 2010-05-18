#ifndef CELL_HH
#define CELL_HH

#include <deque>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "ITrackerGeom/inc/CellId.hh"
#include "ITrackerGeom/inc/CellDetail.hh"
#include "ITrackerGeom/inc/Wire.hh"

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e { 

class Cell{

  friend class SuperLayer;
  friend class ITracker;
  friend class ITrackerMaker;


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

  CLHEP::Hep3Vector getMidPoint() const {return _senseWire.get()->getMidPoint();}

  CLHEP::Hep3Vector getDirection() const { return _senseWire.get()->getDirection();}

//  int hack;
  
protected:

  // Identifier
  CellId _id;

  // Detailed description of a cell.
  boost::shared_ptr<CellDetail> _detail;
  int _detailIndex;

  boost::shared_ptr<Wire> _senseWire;

};

}  //namespace mu2e

#endif /*CELL_HH*/
