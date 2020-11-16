//
// Hold information about one Panel in a tracker.
//
// Original author Rob Kutschke
//

#include <sstream>

#include "TrackerGeom/inc/Panel.hh"

using namespace std;

namespace mu2e {
// define equivalence unique panels
  StrawIdMask Panel::_sidmask(StrawIdMask::uniquepanel);

  string Panel::name( string const& base ) const{
    ostringstream os;
    os << base
       << _id.getPlane() << "_"
       << _id.getPanel();

    return os.str();
  }

  Panel::Panel( const StrawId& id, TrackerStrawCollection const& straws ) : _id(id) {
    for(auto const& straw : straws ) {
    // pick out all the straws belonging to this panel.  This code relies on the Tracker collection being in order.
      if(_sidmask.equal(_id,straw.id())){
	_straws[straw.id().straw()] = &straw;
      }
    }
  }

} // end namespace mu2e
