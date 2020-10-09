
#include "GeneralUtilities/inc/to_string.hh"
#include <iomanip>
#include <sstream>

namespace mu2e{
  std::string to_string( art::SubRunID const& id, int l1, int l2, std::string const& separator){
    std::ostringstream os;
    if ( l1 > 0 && l2 > 0){
      os << std::setfill('0') << std::setw(l1) << id.run() << separator << std::setw(l2) << id.subRun();
    } else{
      os << id.run() << separator << id.subRun();
    }
    return os.str();
  }
}
