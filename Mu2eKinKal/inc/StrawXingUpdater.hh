#ifndef Mu2eKinKal_StrawXingUpdater_hh
#define Mu2eKinKal_StrawXingUpdater_hh
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include <tuple>
namespace mu2e {
  // simple struct to hold crossing calculation configuration parameters
  struct StrawXingUpdater {
    using SXUConfig = std::tuple<float,float,float>;
    WireHitState hitstate_; // associated hit state (inactive if no hit)
    double maxdocasig_; // maximum doca error to use non-averaged value
    double maxdoca_; // maximum DOCA to include this straw's material
    double maxddoca_; // maximum DOCA to use 'exact' calculation, otherwise average over all physical impact parameters
    // default constructor is functional but will always use the impact-parameter averaged material
    StrawXingUpdater() : hitstate_(WireHitState::null), maxdocasig_(-1.0), maxdoca_(0.0), maxddoca_(0.0) {}
    StrawXingUpdater(SXUConfig const& sxusetting) {
      maxdocasig_ = std::get<0>(sxusetting);
      maxdoca_ = std::get<1>(sxusetting);
      maxddoca_ = std::get<2>(sxusetting);
    }
  };
}
#endif
