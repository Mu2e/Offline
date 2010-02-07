//
// Hold information about one Sector in a tracker.
//
//
// $Id: Sector.cc,v 1.1 2010/02/07 00:29:41 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/02/07 00:29:41 $
//
// Original author Rob Kutschke
//

#include <sstream>

#include "TrackerGeom/inc/Sector.hh"

namespace mu2e {

  std::string Sector::name( std::string const& base ) const{
    std::ostringstream os;

    os << base
       << _id.getDevice() << "_"
       << _id.getSector();

    return os.str();
  }

}
