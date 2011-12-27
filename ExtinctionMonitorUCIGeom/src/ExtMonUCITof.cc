//
// Hold information about one Tof (station/segment) in Extinction Monitor.
//
//

#include <sstream>

#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCITof.hh"

#ifndef __CINT__

namespace mu2e {

  namespace ExtMonUCI {

    std::string ExtMonTof::name( std::string const& base ) const {
      std::ostringstream os;

      os << base << "_" << "Tof" << "_"
         << _station << "_" << _segment;
      return os.str();
    }

  } // namespace ExtMonUCI

} // namespace mu2e

#endif

