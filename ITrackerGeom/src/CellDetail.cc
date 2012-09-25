// ITracker cells description
//
// $Id: CellDetail.cc,v 1.7 2012/09/25 10:08:28 tassiell Exp $
// $Author: tassiell $
// $Date: 2012/09/25 10:08:28 $
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
//                        delete *_senseWire; *_senseWire=NULL;
//                } catch (cet::exception e) {
//                    throw cet::exception("GEOM")
//                        << "Error during deleting cell wire detail data \n";
//                }

  }

} // namespace mu2e
