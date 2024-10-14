#ifndef RecoDataProducts_HelixRecoDir_hh
#define RecoDataProducts_HelixRecoDir_hh


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

    // Enum declaration for direction
    enum PropDir {
      upstream = -1,
      ambiguous = 0,
      downstream = 1
    };

    // Method to predict direction based on slope and slopeSig
    PropDir predictDirection(float sigThreshold) const {
      float sig = slopeSig(); // Compute the slope significance

      if (sig < sigThreshold) {
        return ambiguous; // Ambiguous if below the threshold
      } else if (_slope > 0) {
        return downstream; // Downstream if slope > 0
      } else if (_slope < 0) {
        return upstream; // Upstream if slope < 0
      }

      return ambiguous; // Default to ambiguous if slope is 0 or other cases
    }

    //data members
    float _slope;
    float _slopeErr;
    float _chi2ndof;
  };
}
#endif
