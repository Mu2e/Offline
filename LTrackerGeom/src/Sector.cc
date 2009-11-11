//
// Hold information about one Sector in a tracker.
//
//
// $Id: Sector.cc,v 1.1 2009/11/11 14:35:05 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/11 14:35:05 $
//
// Original author Rob Kutschke
//

#include <sstream>

#include "LTrackerGeom/inc/Sector.hh"

namespace mu2e {

  std::string Sector::name( std::string const& base ) const{
    std::ostringstream os;

    os << base
       << _id.getDevice() << "_"
       << _id.getSector();

    return os.str();
  }

}
