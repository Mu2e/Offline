#ifndef ITrackerGeom_ITLayer_hh
#define ITrackerGeom_ITLayer_hh

#include <deque>
#include <vector>

#include <boost/shared_ptr.hpp>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "ITrackerGeom/inc/Cell.hh"
#include "ITrackerGeom/inc/ITLayerId.hh"
#include "ITrackerGeom/inc/ITLayerDetail.hh"

namespace mu2e { 

class ITLayer{

  friend class SuperLayer;
  friend class ITDevice;
  friend class ITracker;
  friend class ITrackerMaker;


public:

  enum Ltype {undefined=-1, wire, gas};

  // A free function, returning void, that takes a const Cell& as an argument.
  typedef void (*CellFunction)( const ITLayer& s);

  ITLayer();

  // Constructor using cell/wire/layer info.
  ITLayer( ITLayerId& id,
         boost::shared_ptr<ITLayerDetail> &detail,
         int detailIndex,
         std::vector<boost::shared_ptr<Cell> > &cells,
         std::vector<boost::shared_ptr<Wire> > &fieldWires
         );
  
  ~ITLayer ();

  const ITLayerId& Id() const { return _id;}

  boost::shared_ptr<ITLayerDetail> getDetail() const { return _detail;}

  const Ltype getLayerType() const { return _layerType; }

  int nCells() const { return _nCells; }

  int nFieldWires() const { return _nFiledWires; }

  double voxelizationFactor() const { return _voxelizationFactor; }

  boost::shared_ptr<Cell> getCell( int n ) const throw(cet::exception) {
        if (n>=0 && n<_nCells) return _cells.at(n);
        else throw cet::exception("GEOM")<< "Cell number: "<< n <<" not present in "<<_id;
  }

  boost::shared_ptr<Cell> getCell( const CellId& id ) const {
    return getCell(id.getCell());
  }

  boost::shared_ptr<Wire> getFWire( int n ) const throw(cet::exception) {
        if (n>=0 && n<_nFiledWires) return _fieldWires.at(n);
        else throw cet::exception("GEOM")<< "Field wire number: "<< n <<" not present in "<<_id;
  }

  boost::shared_ptr<Wire> getFWire( const WireId& id ) const {
    return getFWire(id.getWire());
  }

protected:

  // Identifier
  ITLayerId _id;

  // Detailed description of a layer.
  boost::shared_ptr<ITLayerDetail> _detail;
  int _detailIndex;

  Ltype _layerType;
  int _nCells;
  int _nFiledWires;

  double _voxelizationFactor;

  // Pointers to the cells and field wires in this layer.
  std::vector<boost::shared_ptr<Cell> > _cells;
  std::vector<boost::shared_ptr<Wire> > _fieldWires;

  void addCell(Cell *cell){
          _cells.push_back(boost::shared_ptr<Cell>(cell));
          _nCells++;
  }

  void addFieldWire(Wire *wire){
          _fieldWires.push_back(boost::shared_ptr<Wire> (wire));
          _nFiledWires++;
  }


};

}  //namespace mu2e

#endif /* ITrackerGeom_ITLayer_hh */
