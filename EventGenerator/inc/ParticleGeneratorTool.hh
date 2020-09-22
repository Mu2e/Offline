#ifndef EventGenerator_ParticleGeneratorTool_hh
#define EventGenerator_ParticleGeneratorTool_hh

#include <memory>

#include "art/Utilities/ToolConfigTable.h"

#include "SeedService/inc/SeedService.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"
#include "Mu2eUtilities/inc/GenPhysConfig.hh"

namespace mu2e {

  class ParticleGeneratorTool {
  public:
    virtual void setEngine(art::RandomNumberGenerator::base_engine_t& eng) = 0;
    virtual void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) = 0;
    virtual ~ParticleGeneratorTool() noexcept = default;

    enum SpectrumVar  { TOTAL_ENERGY, KINETIC_ENERY, MOMENTUM };
    static SpectrumVar    parseSpectrumVar(const std::string& name) {
      if (name == "totalEnergy"  )  return TOTAL_ENERGY;
      if (name == "kineticEnergy")  return KINETIC_ENERY;
      if (name == "momentum"     )  return MOMENTUM;
      throw cet::exception("BADCONFIG")<<"ParticleGeneratorTool: unknown spectrum variable "<<name<<"\n";
    }

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Table<GenPhysConfig> physics{Name("physics"), Comment("physics config")};
    };
    typedef art::ToolConfigTable<Config> Parameters;
  };
}

#endif
