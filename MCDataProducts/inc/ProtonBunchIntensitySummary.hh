#ifndef MCDataProducts_ProtonBunchIntensitySummary_hh
#define MCDataProducts_ProtonBunchIntensitySummary_hh
#include <float.h>
#include <math.h>
//
// This object summarizes the proton bunch intensity for a subrun. 
//
// Original author David Brown, LBNL
namespace mu2e {
  class ProtonBunchIntensitySummary {
    public:
      explicit ProtonBunchIntensitySummary() : _nev(0), _iave(0.0), _ivar(0.0) {}
      explicit ProtonBunchIntensitySummary(unsigned nevents, float intensity, float ivar) : _nev(nevents), _iave(intensity), _ivar(ivar) {}
      unsigned count() const { return _nev; }
      float averageIntensity() const { return _iave; }
      float intensityVariance() const { return _ivar; }
      float intensityRMS() const { return sqrt(_ivar); }
      float spillDutyFactor() const { return 1.0/(1.0+_ivar/(_iave*_iave)); }
    private:
      unsigned _nev; // total number of events in this subrun 
      float _iave; // average # of protons/microbunch
      float _ivar; // variance of # of protons/microbunch
  };
}

#endif

