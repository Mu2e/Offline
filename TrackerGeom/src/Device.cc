//
// Hold information about one Device in a tracker.
//
//
// $Id: Device.cc,v 1.1 2010/02/07 00:29:41 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/02/07 00:29:41 $
//
// Original author Rob Kutschke
//

#include <sstream>

#include "TrackerGeom/inc/Device.hh"

namespace mu2e {

  std::string Device::name( std::string const& base ) const{
    std::ostringstream os;

    os << base << _id;

    return os.str();
  }

}
