// Andrei Gaponenko, 2015

#ifndef Mu2eG4_Mu2eG4ResourceLimits_hh
#define Mu2eG4_Mu2eG4ResourceLimits_hh

#include "Mu2eG4/inc/Mu2eG4Config.hh"

namespace mu2e {

  class Mu2eG4ResourceLimits {
    unsigned maxStepsPerTrack_;
    unsigned maxStepPointCollectionSize_;
    unsigned maxSimParticleCollectionSize_;
  public:
    explicit Mu2eG4ResourceLimits(const Mu2eG4Config::Limits& lim);

    unsigned maxStepsPerTrack() const { return maxStepsPerTrack_; }
    unsigned maxStepPointCollectionSize() const { return maxStepPointCollectionSize_; }
    unsigned maxSimParticleCollectionSize() const { return maxSimParticleCollectionSize_; }
  };

} // end namespace mu2e

#endif /* Mu2eG4_Mu2eG4ResourceLimits_hh */
