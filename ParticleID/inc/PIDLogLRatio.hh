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
    LogL hypothesis1_;
    LogL hypothesis2_;
  public:

    double value(double x) const {
      return hypothesis1_.value(x) - hypothesis2_.value(x);
    }

    struct Config {
      fhicl::Table<typename LogL::Config> h1{fhicl::Name("hypothesis1"),
          fhicl::Comment("Distribution for hypothesis 1")
          };

      fhicl::Table<typename LogL::Config> h2 {fhicl::Name("hypothesis2"),
          fhicl::Comment("Distribution for hypothesis 2")
          };
    };

    explicit PIDLogLRatio(const fhicl::ParameterSet& pset);
    explicit PIDLogLRatio(const Config& conf);
  };
}

#endif/*ParticleID_inc_PIDLogLRatio_hh*/
