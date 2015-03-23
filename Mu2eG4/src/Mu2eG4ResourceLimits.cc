#include "Mu2eG4/inc/Mu2eG4ResourceLimits.hh"

#include "fhiclcpp/ParameterSet.h"

namespace mu2e {
  Mu2eG4ResourceLimits::Mu2eG4ResourceLimits(const fhicl::ParameterSet& pset)
    : maxStepsPerTrack_{pset.get<unsigned>("maxStepsPerTrack")}
    , maxStepPointCollectionSize_{pset.get<unsigned>("maxStepPointCollectionSize")}
    , maxSimParticleCollectionSize_{pset.get<unsigned>("maxSimParticleCollectionSize")}
  {}
}
