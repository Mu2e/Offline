// Header describing a Mu2e art::SubRun
//
// Original author Rob Kutschke

#include "Offline/DataProducts/inc/SubRunHeader.hh"
#include <ostream>

namespace mu2e {

  std::ostream& operator<<(std::ostream& os,
			   SubRunHeader const& sh ){
    os << " start Time: " << sh.startTime
       << " first EWT: "  << sh.ewtFirst;
    return os;
  }

}
