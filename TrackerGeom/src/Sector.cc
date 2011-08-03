//
// Hold information about one Sector in a tracker.
//
//
// $Id: Sector.cc,v 1.3 2011/08/03 18:31:25 mf Exp $
// $Author: mf $
// $Date: 2011/08/03 18:31:25 $
//
// Original author Rob Kutschke
//

#include <sstream>

#include "TrackerGeom/inc/Sector.hh"

using namespace std;

namespace mu2e {

  string Sector::name( string const& base ) const{
    ostringstream os;

    os << base
       << _id.getDevice() << "_"
       << _id.getSector();

    return os.str();
  }

  void Sector::fillPointers ( const Tracker& tracker ) const {
    for( size_t i=0; i<_layers.size(); ++i ){
      _layers[i].fillPointers(tracker);
      _straw0MidPoint += _layers[i].straw0MidPoint();
    }
    _straw0Direction = _layers[0].straw0Direction();
    _straw0MidPoint  /= _layers.size();
  }

} // end namespace mu2e
