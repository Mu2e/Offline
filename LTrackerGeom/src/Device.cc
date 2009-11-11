//
// Hold information about one Device in a tracker.
//
//
// $Id: Device.cc,v 1.1 2009/11/11 14:35:05 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/11 14:35:05 $
//
// Original author Rob Kutschke
//

#include <sstream>

#include "LTrackerGeom/inc/Device.hh"

namespace mu2e {

  std::string Device::name( std::string const& base ) const{
    std::ostringstream os;

    os << base << _id;

    return os.str();
  }

}
