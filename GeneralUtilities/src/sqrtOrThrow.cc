//
// Some helper functions for sqrtOrThrow.
//
//
// Contact person Rob Kutschke
//

#include "GeneralUtilities/inc/sqrtOrThrow.hh"

#include <stdexcept>
#include <sstream>

namespace mu2e {

  void throwHelperForSqrtOrThrow( double val, double epsilon){
    std::ostringstream os;
    os << "sqrtOrThrow range error. Inputs " << val << " " << epsilon;
    throw std::range_error(os.str());
  }

  void throwHelperForSqrtOrThrow( float val, float epsilon){
    std::ostringstream os;
    os << "sqrtOrThrow range error. Inputs " << val << " " << epsilon;
    throw std::range_error(os.str());
  }

} // end namespace mu2e
