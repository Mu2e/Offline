#ifndef MCDataProducts_ProtonBunchIntensity_hh
#define MCDataProducts_ProtonBunchIntensity_hh
#include <cmath>
//
// This object defines the proton bunch intensity of a single microbunch (event).  It's
// used when mixing beam-based backgrounds.
//
// Original author David Brown, LBNL
namespace mu2e {
  class ProtonBunchIntensity {
    public:
      explicit ProtonBunchIntensity(unsigned long long intensity=0) : _intensity(intensity) {}
      unsigned long long intensity() const { return _intensity; }
      bool operator == (ProtonBunchIntensity const& other ) const {
        return _intensity == other.intensity(); }
      bool operator != (ProtonBunchIntensity const& other ) const { return !(operator ==(other)); }
      void add(ProtonBunchIntensity const& other) { _intensity += other.intensity(); }
    void set(unsigned long long intensity){ _intensity = intensity;}
    private:
      unsigned long long _intensity; // this has units # of protons/microbunch
  };
}

#endif

