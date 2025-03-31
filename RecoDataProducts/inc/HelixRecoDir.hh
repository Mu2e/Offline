#ifndef RecoDataProducts_HelixRecoDir_hh
#define RecoDataProducts_HelixRecoDir_hh

#include "Offline/RecoDataProducts/inc/TrkFitDirection.hh"

namespace mu2e {

  struct HelixRecoDir {

    HelixRecoDir(float Slope, float SlopeErr, float Chi2ndof):
      _slope(Slope),  _slopeErr(SlopeErr), _chi2ndof(Chi2ndof) {}

    HelixRecoDir():
      _slope(0.),  _slopeErr(0.), _chi2ndof(0.) {}

    //accessors
    float slope()    const { return _slope; }
    float slopeErr() const { return _slopeErr; }
    float slopeSig() const { return std::fabs(_slope/_slopeErr); }
    float chi2ndof() const { return _chi2ndof; }

    // Method to predict direction based on slope and slopeSig
    TrkFitDirection::FitDirection predictDirection(float sigThreshold) const {
      const float sig = slopeSig(); // Compute the slope significance

      if (sig < sigThreshold) {
        return TrkFitDirection::FitDirection::unknown; // Ambiguous if below the threshold
      } else if (_slope > 0.f) {
        return TrkFitDirection::FitDirection::downstream; // Downstream if slope > 0
      } else if (_slope < 0.f) {
        return TrkFitDirection::FitDirection::upstream; // Upstream if slope < 0
      }

      return TrkFitDirection::FitDirection::unknown; // Default to ambiguous if slope is exactly 0 or other cases
    }

    //data members
    float _slope;
    float _slopeErr;
    float _chi2ndof;
  };
}
#endif
