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
    // Physical indicates it followed a physical distribution (including flat if it's physical, e.g. cos(theta_z) for CEs)
    // Flat is for flat and non-physical distributions, like flat electron energy from the target
    enum Type {kPhysical = 0, kFlat, kOther};

    // A variable that can be restricted in a simulation
    struct RestrictedVar {
      RestrictedVar(const std::string name = "default",
                    const double fraction = 1.,
                    const double xmin = std::numeric_limits<double>::lowest(),
                    const double xmax = std::numeric_limits<double>::max(),
                    const Type type = kPhysical) : name_(name),
                                                   fraction_(fraction),
                                                   xmin_(xmin),
                                                   xmax_(xmax),
                                                   type_(type) {}

      // Check if a value is within the restricted range
      bool accepted(double value) const {
        return value >= xmin_ && value <= xmax_; // allow the bounds
      }

      std::string name_     ; // variable identifier, e.g. "energy"
      double      fraction_ ; // normalization correction factor
      double      xmin_     ; // variable range restriction (inclusive)
      double      xmax_     ;
      Type        type_     ; // variable sampling configuration
    };

    SpectrumConfig() {}

    // Add a variable
    void add_var(const RestrictedVar var) {
      vars_[var.name_] = var;
    }

    // Get the variables
    const std::map<std::string, RestrictedVar>& get_variables() const {
      return vars_;
    }

    // Check if provided variables pass the selection
    bool accepted(const std::map<std::string, double>& vars) {
      for(const auto& [name, value] : vars) {
        // check if this name is defined
        if(!vars_.contains(name)) {
          printf("[SpectrumConfig::%s] No variable %s is defined!\n",
                 __func__, name.c_str());
          return false;
        }
        // Check if it's accepted
        if(!vars_[name].accepted(value)) return false;
      }

      // Passes all selections
      return true;
    }

    // Check if all variables are configured as physical
    bool is_physical() const {
      for(const auto& [_, var] : vars_) {
        if(var.type_ != Type::kPhysical) return false;
      }
      return true;
    }

    // Get the total fraction, assuming each component is independent
    double total_fraction(const bool verbose = false) const {
      double fraction = 1.;
      for(const auto& [_,var] : vars_) {
        fraction *= var.fraction_;
        if(verbose) printf("  Component %12s has fraction %.4e\n", var.name_.c_str(), var.fraction_);
      }
      if(verbose) printf("  Total fraction: %.4e\n", fraction);
      return fraction;
    }

    // Get the fraction for one component
    double var_fraction(const std::string name) const {
      // check if this name is defined
      if(!vars_.contains(name)) {
        printf("[SpectrumConfig::%s] No variable %s is defined!\n",
               __func__, name.c_str());
        return 0.;
      }
      return vars_.at(name).fraction_;
    }

  private:
    std::map<std::string,RestrictedVar> vars_; // map of (potentially) restricted variables
  };
}

#endif/*MCDataProducts_inc_SpectrumConfig_hh*/
