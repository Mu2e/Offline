// Relevant event generator configuration information
// Michael MacKenzie, 2026

#ifndef MCDataProducts_inc_SpectrumConfig_hh
#define MCDataProducts_inc_SpectrumConfig_hh

#include <limits>

namespace mu2e {

  class SpectrumConfig {
  public:
    enum Type {kPhysical = 0, kFlat, kOther};

    SpectrumConfig() {}

    // Common operations for normalization
    double ediff()  const { return emax_ - emin_; }
    double czfrac() const { return (czmax_ - czmin_)/2.; }

  public: // allow direct access/manipulation of the fields
    double emin_             = -1.;
    double emax_             = -1.;
    double tmin_             = std::numeric_limits<double>::lowest();
    double tmax_             = std::numeric_limits<double>::max();
    double czmin_            = -1.;
    double czmax_            =  1.;
    double fraction_sampled_ = 1.;
    int    type_             = Type::kOther;
  };
}

#endif/*MCDataProducts_inc_SpectrumConfig_hh*/
