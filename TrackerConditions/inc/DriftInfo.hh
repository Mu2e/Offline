#ifndef Mu2eKinKal_DriftInfo_hh
#define Mu2eKinKal_DriftInfo_hh
#include <algorithm>

namespace mu2e {
  // struct describing drift information
  struct DriftInfo {
    double LorentzAngle_ =0; // angle for EXB effects
    double rDrift_ =0; // calibrated drift distance
    double cDrift_ =0; // single cluster drift distance
    double signedDriftError_ =0; // estimated error on signed drift distance (includes LR ambiguity error effects)
    double unsignedDriftError_ =0; // estimated error on unsigned drift distance
    double driftVelocity_ =0; // instantaneous drift velocity
    static double maxdvar_; // maximum distance variance, given by straw radius
    double driftHitVar() const { return signedDriftError_*signedDriftError_; } // variance for hits constrained to the signed drift distance
    double nullHitVar() const; // variance for hits constrained to the wire position (null hits)
    double unsignedDriftVar() const { return unsignedDriftError_*unsignedDriftError_; }
    double driftTimeVar() const { return unsignedDriftVar()/(driftVelocity_*driftVelocity_); }
  };
}
#endif
