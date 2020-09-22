#ifndef EventGenerator_MuStopParticleGenerator_hh
#define EventGenerator_MuStopParticleGenerator_hh

#include <memory>

#include "art/Utilities/ToolConfigTable.h"

#include "SeedService/inc/SeedService.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"

namespace mu2e {

  class MuStopParticleGenerator {
  public:
    virtual void setEngine(art::RandomNumberGenerator::base_engine_t& eng) = 0;
    virtual void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) = 0;
    virtual ~MuStopParticleGenerator() noexcept = default;

    enum SpectrumVar  { TOTAL_ENERGY, KINETIC_ENERY, MOMENTUM };
    static SpectrumVar    parseSpectrumVar(const std::string& name) {
      if (name == "totalEnergy"  )  return TOTAL_ENERGY;
      if (name == "kineticEnergy")  return KINETIC_ENERY;
      if (name == "momentum"     )  return MOMENTUM;
      throw cet::exception("BADCONFIG")<<"MuStopParticleGenerator: unknown spectrum variable "<<name<<"\n";
    }

    struct Config {
    };
    typedef art::ToolConfigTable<Config> Parameters;
  };
}

#endif
