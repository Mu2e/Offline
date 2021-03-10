#include "Mu2eG4/inc/Mu2eG4Inputs.hh"

namespace mu2e {
  Mu2eG4Inputs::Mu2eG4Inputs(const Mu2eG4Config::Inputs_& conf)
    : primaryType_{conf.primaryType()}
    , primaryTag_{conf.primaryTag()}
    , simParticleNumberOffset_{conf.simParticleNumberOffset()}
    , inputSimParticles_{conf.inputSimParticles()}
    , inputMCTrajectories_{conf.inputMCTrajectories()}
    , inputPhysVolumeMultiInfo_{conf.inputPhysVolumeMultiInfo()}
    , multiStage_{ primaryType_.id() != Mu2eG4PrimaryType::GenParticles }
  {}
}
