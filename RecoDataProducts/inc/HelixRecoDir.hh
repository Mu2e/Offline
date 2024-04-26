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

    //data members
    float _slope;
    float _slopeErr;
    float _chi2ndof;
  };
}
#endif
