// ITracker cells description
//
// $Id: CellDetail.cc,v 1.8 2013/03/15 16:20:00 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 16:20:00 $
//
// Original author G. Tassielli
//

#include "ITrackerGeom/inc/CellDetail.hh"

using namespace std;

namespace mu2e {

  CellDetail::CellDetail( double circumscribedRadius, double inscribedCircleRadius, boost::shared_ptr<WireDetail> senseWire
                           ):
                           _circumscribedRadius(circumscribedRadius),
                           _inscribedCircleRadius(inscribedCircleRadius),
                           _senseWire(senseWire)
  {
  }

  CellDetail::~CellDetail (){
//                try {
//                        delete *_senseWire; *_senseWire=nullptr;
//                } catch (cet::exception e) {
//                    throw cet::exception("GEOM")
//                        << "Error during deleting cell wire detail data \n";
//                }

  }

} // namespace mu2e
