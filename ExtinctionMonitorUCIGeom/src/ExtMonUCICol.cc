//
// Hold information about one Collimator in Extinction Monitor.
//
//
// $Id: ExtMonUCICol.cc,v 1.1 2011/12/10 00:16:15 youzy Exp $
// $Author: youzy $
// $Date: 2011/12/10 00:16:15 $
//

#include <sstream>

#include "ExtinctionMonitorUCIGeom/inc/ExtMonUCICol.hh"

#ifndef __CINT__

namespace mu2e {

  namespace ExtMonUCI {

    std::string ExtMonCol::name( std::string const& base ) const {
      std::ostringstream os;

      os << base << "_" << "Col" << "_"
         << _id;
      return os.str();
    }

  } // namespace ExtMonUCI

} // namespace mu2e

#endif

