// Framework includes
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// Mu2e includes
#include "ITrackerGeom/inc/Cell.hh"

#ifndef __CINT__ 

using CLHEP::Hep3Vector;

namespace mu2e {

Cell::Cell():
  _id(CellId()),
  _detailIndex(0)
{
}

Cell::Cell( CellId id,
		 boost::shared_ptr<CellDetail> detail,
		 boost::shared_ptr<Wire> senseWire
		):
  _id(id),
  _detail(detail),
  _senseWire(senseWire)
{
}

Cell::~Cell (){
//	try {
//		if (_detail) delete _detail;
//		if (_senseWire) delete _senseWire;
//	} catch (cms::Exception e) {
//	    throw cms::Exception("GEOM")
//		<< "Error during deleting cell data \n";
//	}
}

} // namespace mu2e

#endif
