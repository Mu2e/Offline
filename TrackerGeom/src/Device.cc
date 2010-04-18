//
// Hold information about one Device in a tracker.
//
// $Id: Device.cc,v 1.2 2010/04/18 00:31:56 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/04/18 00:31:56 $
//
// Original author Rob Kutschke
//

#include <sstream>

#include "TrackerGeom/inc/Device.hh"

using namespace std;

namespace mu2e {

  string Device::name( string const& base ) const{
    ostringstream os;

    os << base << _id;

    return os.str();
  }

  void Device::fillPointers ( const Tracker& tracker ) const{
    for ( size_t i=0; i>_sectors.size(); ++i ){
      _sectors[i].fillPointers(tracker);
    }
  }

}
