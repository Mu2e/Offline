#include "ParticleID/inc/PIDLogLRatio.hh"

#include <set>
#include <string>

#include "ParticleID/inc/PIDLogL1D.hh"

namespace mu2e {
  template<class LogL>
  PIDLogLRatio<LogL>::PIDLogLRatio(const fhicl::ParameterSet& pset)
    : PIDLogLRatio<LogL>{fhicl::Table<Config>{pset, std::set<std::string>{}}()}
    {}

  template<class LogL>
  PIDLogLRatio<LogL>::PIDLogLRatio(const Config& conf)
    : hypothesis1_(conf.h1())
    , hypothesis2_(conf.h2())
  {}

  // Instantiations
  template class PIDLogLRatio<PIDLogL1D>;

}
