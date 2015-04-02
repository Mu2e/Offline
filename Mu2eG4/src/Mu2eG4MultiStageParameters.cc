#include "Mu2eG4/inc/Mu2eG4MultiStageParameters.hh"

#include "fhiclcpp/ParameterSet.h"

namespace mu2e {
  Mu2eG4MultiStageParameters::Mu2eG4MultiStageParameters(const fhicl::ParameterSet& pset)
    : multiStage_(!pset.is_empty())
    , simParticleNumberOffset_{multiStage_ ? pset.get<unsigned>("simParticleNumberOffset") : 0 }
    , inputSimParticles_{multiStage_ ? pset.get<std::string>("inputSimParticles") : "" }
    , inputMCTrajectories_{multiStage_ ? pset.get<std::string>("inputMCTrajectories") : "" }
    , inputPhysVolumeMultiInfo_{multiStage_ ? pset.get<std::string>("inputPhysVolumeMultiInfo") : "" }
    , genInputHits_{multiStage_ ? pset.get<std::vector<art::InputTag> >("genInputHits") : std::vector<art::InputTag>() }
  {}
}
