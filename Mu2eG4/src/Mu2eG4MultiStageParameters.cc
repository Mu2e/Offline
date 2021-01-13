#include "Mu2eG4/inc/Mu2eG4MultiStageParameters.hh"

namespace {
  // Because a Mu2eG4MultiStageParameters object can be created from
  // multiple threads and because it uses configuration validation,
  // which is not thread-safe, there exists the possibility of data
  // races.
  //
  // To avoid this, we create a global configuration struct, which
  // does all of the thread-unsafe manipulations upon construction
  // (namely adjusting FHiCL's table-member and name-stack
  // registries). We then copy the struct in the body of the c'tor,
  // where they copy operation does not use the name-stack registry.
  //
  // NB - This is an expert-only workaround which should arguably go
  //      away if fhiclcpp decides to adopt thread-local registries
  //      for configuration validation.
  mu2e::Mu2eG4Config::MultiStageParameters_ const global_msp;
}

namespace mu2e {
  Mu2eG4MultiStageParameters::Mu2eG4MultiStageParameters(const Mu2eG4Config::Top& conf)
    : multiStage_(false)
    , simParticleNumberOffset_{0 }
    , inputSimParticles_{""}
    , inputMCTrajectories_{""}
    , inputPhysVolumeMultiInfo_{""}
    , genInputHits_{std::vector<art::InputTag>()}
  {
    auto msp = global_msp;
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
