#include "Mu2eG4/inc/Mu2eG4ResourceLimits.hh"

namespace mu2e {
  Mu2eG4ResourceLimits::Mu2eG4ResourceLimits(const Mu2eG4Config::Limits& lim)
    : maxStepsPerTrack_{lim.maxStepsPerTrack()}
    , maxStepPointCollectionSize_{lim.maxStepPointCollectionSize()}
    , maxSimParticleCollectionSize_{lim.maxSimParticleCollectionSize()}
  {}
}
