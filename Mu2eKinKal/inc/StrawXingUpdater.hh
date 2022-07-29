#ifndef Mu2eKinKal_StrawXingUpdater_hh
#define Mu2eKinKal_StrawXingUpdater_hh
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include <tuple>
namespace mu2e {
  // simple struct to hold crossing calculation configuration parameters
  struct StrawXingUpdater {
    using SXUConfig = std::tuple<float,float,float,float,bool>;
    WireHitState hitstate_; // associated hit state (inactive if no hit)
    double maxdoca_; // maximum DOCA to include this straw's material
    double maxddoca_; // maximum DOCA to use non-averaged value
    double maxdocasig_; // maximum doca error to use non-averaged value
    double varscale_; // scaling factor for variance
    bool scalevar_; // scale variance or not
    // default constructor is functional but will always use the impact-parameter averaged material
    StrawXingUpdater() : hitstate_(WireHitState::null), maxdoca_(0.0), maxddoca_(0.0), maxdocasig_(-1.0),  varscale_(1.0), scalevar_(true) {}
    StrawXingUpdater(SXUConfig const& sxusetting) {
      maxdoca_ = std::get<0>(sxusetting);
      maxddoca_ = std::get<1>(sxusetting);
      maxdocasig_ = std::get<2>(sxusetting);
      varscale_ = std::get<3>(sxusetting);
      scalevar_ = std::get<4>(sxusetting);
    }
  };
}
#endif
