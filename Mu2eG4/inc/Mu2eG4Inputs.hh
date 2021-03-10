// Andrei Gaponenko, 2015, 2021

#ifndef Mu2eG4_Mu2eG4Inputs_hh
#define Mu2eG4_Mu2eG4Inputs_hh

#include <string>
#include <vector>

#include "canvas/Utilities/InputTag.h"
#include "Mu2eG4/inc/Mu2eG4Config.hh"
#include "Mu2eG4/inc/Mu2eG4PrimaryType.hh"

namespace mu2e {

  class Mu2eG4Inputs {
  public:
    explicit Mu2eG4Inputs(const Mu2eG4Config::Inputs_& conf);

    bool multiStage() const { return multiStage_; }

    Mu2eG4PrimaryType primaryType() const { return primaryType_; }

    const art::InputTag& primaryTag() const { return primaryTag_; }

    unsigned simParticleNumberOffset() const { return simParticleNumberOffset_; }

    const art::InputTag& inputSimParticles() const { return inputSimParticles_; }

    const art::InputTag& inputMCTrajectories() const { return inputMCTrajectories_; }

    const art::InputTag& inputPhysVolumeMultiInfo() const { return inputPhysVolumeMultiInfo_; }

  private:
    Mu2eG4PrimaryType primaryType_;
    art::InputTag primaryTag_;
    unsigned simParticleNumberOffset_;
    art::InputTag inputSimParticles_;
    art::InputTag inputMCTrajectories_;
    art::InputTag inputPhysVolumeMultiInfo_;

    // derived
    bool multiStage_;
  };

} // end namespace mu2e

#endif /* Mu2eG4_Mu2eG4Inputs_hh */
