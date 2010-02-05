#ifndef CELLDETAIL_HH
#define CELLDETAIL_HH

#include <vector>
#include <string>

#include <boost/shared_ptr.hpp>

// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"

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
  double      wireRadius() const throw(cms::Exception) {
	  try {
		  return _senseWire.get()->outerRadius();
	  } catch (cms::Exception e) {
	      throw cms::Exception("GEOM")
		<< "No sense wire defined for the Cell \n";
          return 0.0;
	  }
  }

  double      halfLength() const throw(cms::Exception) {
	  try {
		  return _senseWire->halfLength();
	  } catch (cms::Exception e) {
	      throw cms::Exception("GEOM")
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

#endif /*CELLDETAIL_HH*/
