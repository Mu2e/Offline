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
//                        if (*_senseWire) {delete *_senseWire; *_senseWire=NULL;};
//                } catch (cet::exception e) {
//                    throw cet::exception("GEOM")
//                        << "Error during deleting cell wire detail data \n";
//                }

  }

} // namespace mu2e
