//
// Hold information about one Shielding block in Extinction Monitor.
//
// $Id: ExtMonUCIShd.cc,v 1.1 2012/02/16 20:25:46 youzy Exp $
// $Author: youzy $
// $Date: 2012/02/16 20:25:46 $

#include <sstream>

#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCIShd.hh"

#ifndef __CINT__

namespace mu2e {

  namespace ExtMonUCI {

    std::string ExtMonShd::name( std::string const& base ) const {
      std::ostringstream os;

      os << base << "_" << "Shd" << "_"
         << _id;
      return os.str();
    }

  } // namespace ExtMonUCI

} // namespace mu2e

#endif

