// Header describing a Mu2e art::Event
//
// Original author Rob Kutschke

#include "Offline/DataProducts/inc/EventHeader.hh"
#include <ostream>

namespace mu2e {

  std::ostream& operator<<(std::ostream& os,
			   EventHeader const& eh ){
    os <<  "EWT: "           << eh.ewt
       << " Mode: "          << eh.mode
       << " RF Marker TDC: " << eh.rfmTDC
       << " Flags: "         << eh.flags;
    return os;
  }

}
