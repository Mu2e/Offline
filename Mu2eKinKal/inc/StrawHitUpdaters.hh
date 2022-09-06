//
// enum of StrawHitUpdaters.  This needs to be extended when new updaters are added
//
#ifndef Mu2eKinKal_StrawHitUpdaters_hh
#define Mu2eKinKal_StrawHitUpdaters_hh
#include <string>
#include <vector>

namespace mu2e {
  // straw hit updater algorithms: this needs to be extended if new updaters are defined
  struct StrawHitUpdaters {
    enum algorithm: int {none=0, CA=1, ANN=2, Bkg=3, Combinatoric=10, nalgos };
    // translate from algo to name
    static std::string const& name(algorithm alg);
    static algorithm algo(std::string const& name);
    static std::vector<std::string> names_;
  };
}
#endif
