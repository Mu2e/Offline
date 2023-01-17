#ifndef Mu2eKinKal_DriftInfo_hh
#define Mu2eKinKal_DriftInfo_hh
#include <algorithm>

namespace mu2e {
  // struct describing drift information
  struct DriftInfo {
    double LorentzAngle_ =0; // angle for EXB effects
    double rDrift_ =0; // calibrated drift distance
    double cDrift_ =0; // single cluster drift distance
    double signedDriftError_ =0; // estimated variance on drift distance (includes LR ambiguity errors)
    double unsignedDriftError_ =0; // estimated variance on drift distance for null hits (unsigned)
    double driftVelocity_; // instantaneous drift velocity
    static double maxdvar_; // maximum distance variance, given by straw radius
    double signedDriftVar() const { return signedDriftError_*signedDriftError_; }
    double unsignedDriftVar() const;
  };
}
#endif
