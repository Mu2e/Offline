// Andrei Gaponenko, 2015

#ifndef Mu2eG4_Mu2eG4MultiStageParameters_hh
#define Mu2eG4_Mu2eG4MultiStageParameters_hh

#include <string>
#include <vector>

#include "canvas/Utilities/InputTag.h"

#include "Mu2eG4/inc/Mu2eG4Config.hh"

namespace mu2e {

  class Mu2eG4MultiStageParameters {
  public:
    explicit Mu2eG4MultiStageParameters(const Mu2eG4Config::Top& conf);

    bool multiStage() const { return multiStage_; }

    unsigned simParticleNumberOffset() const { return simParticleNumberOffset_; }

    const art::InputTag& inputSimParticles() const { return inputSimParticles_; }

    const art::InputTag& inputMCTrajectories() const { return inputMCTrajectories_; }

    const art::InputTag& inputPhysVolumeMultiInfo() const { return inputPhysVolumeMultiInfo_; }

    const std::vector<art::InputTag>& genInputHits() const { return genInputHits_; }

  private:
    bool multiStage_;
    unsigned simParticleNumberOffset_;
    art::InputTag inputSimParticles_;
    art::InputTag inputMCTrajectories_;
    art::InputTag inputPhysVolumeMultiInfo_;
    std::vector<art::InputTag> genInputHits_;
  };

} // end namespace mu2e

#endif /* Mu2eG4_Mu2eG4MultiStageParameters_hh */
