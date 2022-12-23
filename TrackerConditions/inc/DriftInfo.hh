#ifndef Mu2eKinKal_DriftInfo_hh
#define Mu2eKinKal_DriftInfo_hh
#include <algorithm>

namespace mu2e {
  // struct describing drift information
  struct DriftInfo {
    double LorentzAngle_; // angle for EXB effects
    double driftDistance_; // drift distance, calculated from drift time
    double driftDistanceError_; // estimated variance on drift distance
    double driftVelocity_; // instantaneous drift velocity
    static double maxdvar_; // maximum distance variance, given by straw radius
    double driftDistanceVar() const { return driftDistanceError_*driftDistanceError_; }
    double nullDistanceVar() const {
      return std::min(maxdvar_, driftDistance_*driftDistance_ + driftDistanceVar()); }
  };
}
#endif
