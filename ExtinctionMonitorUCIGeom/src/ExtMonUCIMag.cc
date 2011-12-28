//
// Hold information about one Magnet in Extinction Monitor.
//
// $Id: ExtMonUCIMag.cc,v 1.2 2011/12/28 00:25:05 youzy Exp $
// $Author: youzy $
// $Date: 2011/12/28 00:25:05 $

#include <sstream>

#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCIMag.hh"

#ifndef __CINT__

namespace mu2e {

  namespace ExtMonUCI {

    std::string ExtMonMag::name( std::string const& base ) const {
      std::ostringstream os;

      os << base << "_" << "Mag" << "_"
         << _id;
      return os.str();
    }

  } // namespace ExtMonUCI

} // namespace mu2e

#endif

