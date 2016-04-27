#include "ParticleID/inc/PIDLogLRatio.hh"

#include <set>
#include <string>

#include "ParticleID/inc/PIDLogL1D.hh"
#include "ParticleID/inc/PIDLogLEp.hh"

namespace mu2e {
  template<class LogL>
  PIDLogLRatio<LogL>::PIDLogLRatio(const fhicl::ParameterSet& pset)
    : PIDLogLRatio<LogL>{fhicl::Table<Config>{pset, std::set<std::string>{}}()}
    {}

  template<class LogL>
  PIDLogLRatio<LogL>::PIDLogLRatio(const Config& conf)
    : signalHypothesis_(conf.signalHypothesis())
    , backgroundHypothesis_(conf.backgroundHypothesis())
  {}

  // Instantiations
  template class PIDLogLRatio<PIDLogL1D>;
  template class PIDLogLRatio<PIDLogLEp>;
}
