#ifndef MCDataProducts_ProtonBunchIntensity_hh
#define MCDataProducts_ProtonBunchIntensity_hh
#include <float.h>
#include <math.h>
//
// This object defines the proton bunch intensity of a single microbunch (event).  It's
// used when mixing beam-based backgrounds.
//
// Original author David Brown, LBNL
namespace mu2e {
  class ProtonBunchIntensity {
    public:
      explicit ProtonBunchIntensity(unsigned intensity=0) : _intensity(intensity) {}
      unsigned intensity() const { return _intensity; }
      bool operator == (ProtonBunchIntensity const& other ) const {
	return _intensity == other.intensity(); }
      bool operator != (ProtonBunchIntensity const& other ) const { return !(operator ==(other)); }
      void add(ProtonBunchIntensity const& other) { _intensity += other.intensity(); }
    private:
      unsigned _intensity; // this has units # of protons/microbunch
  };
}

#endif

