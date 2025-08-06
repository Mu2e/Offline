#ifndef RecoDataProducts_RecoProtonBunchIntensity_hh
#define RecoDataProducts_RecoProtonBunchIntensity_hh
#include <cmath>
//
// This object defines the reconstructed proton bunch intensity of a single microbunch (event).
//
// Original author M. MacKenzie, 2025
namespace mu2e {
  class RecoProtonBunchIntensity {
    public:
    explicit RecoProtonBunchIntensity(unsigned long long intensity = 0, unsigned long long uncertainty = 0) :
      _intensity(intensity), _uncertainty(uncertainty) {}
      unsigned long long intensity() const { return _intensity; }
      unsigned long long uncertainty() const { return _uncertainty; }
      bool operator == (RecoProtonBunchIntensity const& other ) const {
        return _intensity == other.intensity() && _uncertainty == other.uncertainty(); }
      bool operator != (RecoProtonBunchIntensity const& other ) const { return !(operator ==(other)); }
    void add(RecoProtonBunchIntensity const& other) { _intensity += other.intensity(); _uncertainty = std::sqrt(std::pow(_uncertainty, 2) + std::pow(other.uncertainty(), 2)); }
    private:
      unsigned long long _intensity; // this has units # of protons/microbunch
      unsigned long long _uncertainty;
  };
}

#endif
