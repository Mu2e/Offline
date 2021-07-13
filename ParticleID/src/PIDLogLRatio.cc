#include "Offline/ParticleID/inc/PIDLogLRatio.hh"

#include <set>
#include <string>

#include "Offline/ParticleID/inc/PIDLogL1D.hh"
#include "Offline/ParticleID/inc/PIDLogLEp.hh"

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
