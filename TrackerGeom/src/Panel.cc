//
// Hold information about one Panel in a tracker.
//
//
// $Id: Panel.cc,v 1.3 2011/08/03 18:31:25 mf Exp $
// $Author: mf $
// $Date: 2011/08/03 18:31:25 $
//
// Original author Rob Kutschke
//

#include <sstream>

#include "TrackerGeom/inc/Panel.hh"

using namespace std;

namespace mu2e {

  string Panel::name( string const& base ) const{
    ostringstream os;

    os << base
       << _id.getPlane() << "_"
       << _id.getPanel();

    return os.str();
  }

  void Panel::fillPointers ( const Tracker& tracker ) const {
    for( size_t i=0; i<_layers.size(); ++i ){
      _layers[i].fillPointers(tracker);
      _straw0MidPoint += _layers[i].straw0MidPoint();
    }
    _straw0Direction = _layers[0].straw0Direction();
    _straw0MidPoint  /= _layers.size();
  }

} // end namespace mu2e
