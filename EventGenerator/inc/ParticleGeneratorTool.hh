#ifndef EventGenerator_ParticleGeneratorTool_hh
#define EventGenerator_ParticleGeneratorTool_hh

#include <memory>
#include <vector>

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

    virtual void setEngine(art::RandomNumberGenerator::base_engine_t& eng) = 0;
    virtual void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) = 0;
    virtual std::vector<Kinematic> generate() = 0;

    virtual ~ParticleGeneratorTool() noexcept = default;

    enum SpectrumVar  { TOTAL_ENERGY, KINETIC_ENERY, MOMENTUM };
    static SpectrumVar    parseSpectrumVar(const std::string& name) {
      if (name == "totalEnergy"  )  return TOTAL_ENERGY;
      if (name == "kineticEnergy")  return KINETIC_ENERY;
      if (name == "momentum"     )  return MOMENTUM;
      throw cet::exception("BADCONFIG")<<"ParticleGeneratorTool: unknown spectrum variable "<<name<<"\n";
    }
  };
}

#endif
