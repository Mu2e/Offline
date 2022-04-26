// Header describing a Mu2e art::Run
//
// Original author Rob Kutschke

#include "Offline/DataProducts/inc/RunHeader.hh"
#include <ostream>

namespace mu2e {

  std::ostream& operator<<(std::ostream& os,
			   RunHeader const& rh ){
    os << " start Time: " << rh.startTime
       << " first EWT: "  << rh.ewtFirst;
    return os;
  }

}
