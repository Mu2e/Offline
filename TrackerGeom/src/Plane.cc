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

  string Plane::name( string const& base ) const{
    ostringstream os;

    os << base << _id;

    return os.str();
  }

  void Plane::fillPointers ( const Tracker* tracker ) const{
    for ( size_t i=0; i<_panels.size(); ++i ){
      _panels[i]->fillPointers(tracker);
    }
  }

}
