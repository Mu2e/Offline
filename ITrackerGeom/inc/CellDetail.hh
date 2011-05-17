#ifndef ITrackerGeom_CellDetail_hh
#define ITrackerGeom_CellDetail_hh

#include <vector>
#include <string>

#include <boost/shared_ptr.hpp>

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "ITrackerGeom/inc/WireDetail.hh"

namespace mu2e {

class CellDetail{

public:
  CellDetail():
          _circumscribedRadius(0.0),
          _inscribedCircleRadius(0.0)
  {}

  CellDetail( double circumscribedRadius, double inscribedCircleRadius, boost::shared_ptr<WireDetail> senseWire );
  
  ~CellDetail ();

  double      CirumscribedRadius()   const { return _circumscribedRadius;}
  double      InscribedCircleRadius()   const { return _inscribedCircleRadius;}
  double      wireRadius() const throw(cet::exception) {
          try {
                  return _senseWire.get()->outerRadius();
          } catch (cet::exception e) {
              throw cet::exception("GEOM")
                << "No sense wire defined for the Cell \n";
          return 0.0;
          }
  }

  double      halfLength() const throw(cet::exception) {
          try {
                  return _senseWire->halfLength();
          } catch (cet::exception e) {
              throw cet::exception("GEOM")
                << "No sense wire defined for the Cell \n";
          return 0.0;
          }

  }

private:

  double _circumscribedRadius;
  double _inscribedCircleRadius;
  boost::shared_ptr<WireDetail> _senseWire;

};

}  //namespace mu2e

#endif /* ITrackerGeom_CellDetail_hh */
