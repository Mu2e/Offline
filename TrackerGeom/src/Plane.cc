//
// Hold information about one Plane in a tracker.
//
//
// Original author Rob Kutschke
//

#include <sstream>

#include "TrackerGeom/inc/Plane.hh"

using namespace std;

namespace mu2e {
// define equivalence plane 
  StrawIdMask Plane::_sidmask(StrawIdMask::plane);

  string Plane::name( string const& base ) const{
    ostringstream os;

    os << base << _id;

    return os.str();
  }

  Plane::Plane( const StrawId& id, TrackerPanelCollection const& panels) : _id(id) {
    for(auto const& panel : panels ) {
    // pick out all the panels belonging to this plane.  This code relies on the Tracker collection being in order.
      if(_sidmask.equal(_id,panel.id())){
	_panels[panel.id().panel()] = &panel;
      }
    }
  }

}
