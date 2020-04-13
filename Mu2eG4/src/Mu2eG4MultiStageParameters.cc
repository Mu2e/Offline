#include "Mu2eG4/inc/Mu2eG4MultiStageParameters.hh"

namespace mu2e {
  Mu2eG4MultiStageParameters::Mu2eG4MultiStageParameters(const Mu2eG4Config::Top& conf)
    : multiStage_(false)
    , simParticleNumberOffset_{0 }
    , inputSimParticles_{""}
    , inputMCTrajectories_{""}
    , inputPhysVolumeMultiInfo_{""}
    , genInputHits_{std::vector<art::InputTag>()}
  {
    Mu2eG4Config::MultiStageParameters_ msp;
    multiStage_ = conf.MultiStageParameters(msp);
    if(multiStage_) {
      simParticleNumberOffset_ = msp.simParticleNumberOffset();
      inputSimParticles_ = msp.inputSimParticles();
      inputMCTrajectories_ = msp.inputMCTrajectories();
      inputPhysVolumeMultiInfo_ = msp.inputPhysVolumeMultiInfo();
      genInputHits_ = msp.genInputHits();
    }

  }
}
