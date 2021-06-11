#ifndef EventGenerator_ParticleGeneratorTool_hh
#define EventGenerator_ParticleGeneratorTool_hh

#include <memory>
#include <vector>
#include <string>

#include "CLHEP/Vector/LorentzVector.h"

#include "art/Utilities/ToolConfigTable.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/ProcessCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"

#include "GeneralUtilities/inc/RSNTIO.hh"

namespace mu2e {

  class ParticleGeneratorTool {
  public:

    struct Kinematic {
      PDGCode::type pdgId;
      ProcessCode creationCode;
      CLHEP::HepLorentzVector fourmom;
    };

    virtual void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string& materialName) = 0;

    virtual std::vector<Kinematic> generate() = 0;

    // This interface should be removed when we retire ntuple-based muon resampling
    virtual void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) = 0;

    virtual ~ParticleGeneratorTool() noexcept = default;
  };
}

#endif
