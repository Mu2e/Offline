//
// Hold information about one Tof (station/segment) in Extinction Monitor.
//
// $Id: ExtMonUCITof.cc,v 1.2 2011/12/28 00:25:05 youzy Exp $
// $Author: youzy $
// $Date: 2011/12/28 00:25:05 $

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

