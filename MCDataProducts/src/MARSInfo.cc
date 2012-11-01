#include "MCDataProducts/inc/MARSInfo.hh"

namespace mu2e {
  std::ostream& operator<<(std::ostream& os, const MARSInfo& mi) {
    return os<<"MARSInfo(proton="<<mi.protonNumber()
             <<", subRun="<<mi.subRunNumber()
             <<", run="<<mi.runNumber()
             <<", weight="<<mi.weight()
             <<" )";
  }
}
