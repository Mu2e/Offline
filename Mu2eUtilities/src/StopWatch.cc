#include "Offline/Mu2eUtilities/inc/StopWatch.hh"

namespace mu2e {
  std::ostream& operator<<(std::ostream& os, const StopWatch& watch) {
    watch.Print(os);
    return os;
  }
}
