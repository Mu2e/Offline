#include "MCDataProducts/inc/TrackSummaryTruthAssns.hh"
#include <ostream>

namespace mu2e {
  std::ostream& operator<<(std::ostream& os, const TrackSummaryMatchInfo& mi) {
    return os<<"TrackSummaryMatchInfo(np="<<mi.nPrincipal()
             <<", nAll="<<mi.nAll()
             <<" )";
  }
}
