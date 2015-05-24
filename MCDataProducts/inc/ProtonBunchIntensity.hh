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
      ProtonBunchIntensity(unsigned intensity=0,double meanIntensity=0) : _intensity(intensity),
	  _meanIntensity(meanIntensity){}
      unsigned intensity() const { return _intensity; }
      float meanIntensity() const { return _meanIntensity; }
      double weight() const { return _intensity/_meanIntensity; }
      bool operator == (ProtonBunchIntensity const& other ) const {
	return _intensity == other.intensity() && 
	  fabs(_meanIntensity-other.meanIntensity()) < FLT_MIN; }
      bool operator != (ProtonBunchIntensity const& other ) const { return !(operator ==(other)); }
      void add(ProtonBunchIntensity const& other) { _intensity += other.intensity();
						    _meanIntensity += other.meanIntensity(); }
    private:
      unsigned _intensity; // this has units # of protons/microbunch
      float _meanIntensity; // mean intensity used to generate this bunch
  };
}

#endif

