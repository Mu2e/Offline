#ifndef RecoDataProducts_RecoProtonBunchIntensity_hh
#define RecoDataProducts_RecoProtonBunchIntensity_hh
//
// This object defines the reconstructed proton bunch intensity of a single microbunch (event).
//
// Original author M. MacKenzie, 2025

#include <cmath>
#include <functional>

namespace mu2e {
  class RecoProtonBunchIntensity {
    public:
    typedef unsigned long long RecoPBI_t;
    explicit RecoProtonBunchIntensity(RecoPBI_t intensity = 0, RecoPBI_t uncertainty = 0) :
      _intensity(intensity), _uncertainty(uncertainty) {}
      RecoPBI_t intensity() const { return _intensity; }
      RecoPBI_t uncertainty() const { return _uncertainty; }
      bool operator == (RecoProtonBunchIntensity const& other ) const {
        return _intensity == other.intensity() && _uncertainty == other.uncertainty(); }
      bool operator != (RecoProtonBunchIntensity const& other ) const { return !(operator ==(other)); }
    void setIntensity(RecoPBI_t tmp) { _intensity = tmp; }
    void setUncertainty(RecoPBI_t tmp) { _uncertainty = tmp; }
    template<typename UncAdd> void add(RecoProtonBunchIntensity const& other, UncAdd unc_add = unc_add_quad) {
      _intensity += other.intensity();
      _uncertainty = unc_add(_uncertainty, other.uncertainty());
    }

    // common uncertainty combination options
    static RecoPBI_t unc_add_quad(const RecoPBI_t u_1, const RecoPBI_t u_2) { return std::sqrt(std::pow(u_1, 2) + std::pow(u_2, 2)); }
    static RecoPBI_t unc_add_linear(const RecoPBI_t u_1, const RecoPBI_t u_2) { return u_1 + u_2; }
    private:
      RecoPBI_t _intensity; // this has units # of protons/microbunch
      RecoPBI_t _uncertainty;
  };
}

#endif
