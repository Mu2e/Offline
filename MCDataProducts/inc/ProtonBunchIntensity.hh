#ifndef MCDataProducts_ProtonBunchIntensity_hh
#define MCDataProducts_ProtonBunchIntensity_hh
//
// This object defines the proton bunch intensity of a single microbunch (event).  It's
// used when mixing beam-based backgrounds.
//
// Original author David Brown, LBNL
namespace mu2e {
  class ProtonBunchIntensity {
    public:
      ProtonBunchIntensity(unsigned intensity=0) : _intensity(intensity) {}
      unsigned intensity() const { return _intensity; }
    private:
      unsigned _intensity; // this has units # of protons/microbunch
  };
}

#endif

