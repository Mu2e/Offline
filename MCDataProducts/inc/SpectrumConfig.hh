// Relevant event generator configuration information
// Michael MacKenzie, 2026

#ifndef MCDataProducts_inc_SpectrumConfig_hh
#define MCDataProducts_inc_SpectrumConfig_hh

#include <limits>
#include <vector>
#include <string>

namespace mu2e {

  class SpectrumConfig {
  public:
    // Indicate the broad class of primary simulations this falls under
    enum Type {kPhysical = 0, kFlat, kOther};

    // A variable that can be restricted in a simulation
    struct RestrictedVar {
      RestrictedVar(const std::string name = "default",
                    const double fraction = 1.,
                    const double xmin = std::numeric_limits<double>::lowest(),
                    const double xmax = std::numeric_limits<double>::max()) : name_(name),
                                                                              fraction_(fraction),
                                                                              xmin_(xmin),
                                                                              xmax_(xmax) {}

      std::string name_     ;
      double      fraction_ ;
      double      xmin_     ;
      double      xmax_     ;
    };

    SpectrumConfig(const int type = Type::kOther) : type_(type) {}

    double total_fraction(const bool verbose = false) {
      double fraction = 1.;
      for(const auto& var : vars_) {
        fraction *= var.fraction_;
        if(verbose) printf("  Component %12s has fraction %.4e\n", var.name_.c_str(), var.fraction_);
      }
      if(verbose) printf("  Total fraction: %.4e\n", fraction);
      return fraction;
    }

  public: // allow direct access/manipulation of the fields
    std::vector<RestrictedVar> vars_;
    int type_;
  };
}

#endif/*MCDataProducts_inc_SpectrumConfig_hh*/
