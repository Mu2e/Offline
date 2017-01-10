// A log likelihood ratio.
//
// Andrei Gaponenko, 2016

#ifndef ParticleID_inc_PIDLogLRatio_hh
#define ParticleID_inc_PIDLogLRatio_hh

#include "fhiclcpp/types/Table.h"

#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

  template<class LogL>
  class PIDLogLRatio {
    LogL signalHypothesis_;
    LogL backgroundHypothesis_;
  public:

    template<typename ... Args>
    double value(Args&&... x) const {
      const double ps = signalHypothesis_.value(x...);
      // Return cutoff() if the inputs are inconsistent with
      // signal, no matter what their consistency with background is.
      return (ps <= LogL::cutoff()) ?
        LogL::cutoff() :
        ps - backgroundHypothesis_.value(x...);
    }

    static double cutoff() { return LogL::cutoff(); }

    const LogL& signalLogLikelihood() const { return signalHypothesis_; }
    const LogL& backgroundLogLikelihood() const { return backgroundHypothesis_; }

    struct Config {
      fhicl::Table<typename LogL::Config> signalHypothesis{fhicl::Name("signalHypothesis"),
          fhicl::Comment("Distribution for signal hypothesis")
          };

      fhicl::Table<typename LogL::Config> backgroundHypothesis{fhicl::Name("backgroundHypothesis"),
          fhicl::Comment("Distribution for background hypothesis")
          };
    };

    explicit PIDLogLRatio(const fhicl::ParameterSet& pset);
    explicit PIDLogLRatio(const Config& conf);
  };
}

#endif/*ParticleID_inc_PIDLogLRatio_hh*/
