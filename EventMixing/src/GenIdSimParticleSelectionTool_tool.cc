// Ed Callaghan
// Tool to select SimParticles downstread by a specific Mu2e generator
// November 2024

#include "Offline/EventMixing/inc/GenIdSimParticleSelectionTool.hh"

namespace mu2e{

  GenIdSimParticleSelectionTool::GenIdSimParticleSelectionTool(const Parameters& config):
      _genId(GenId::findByName(config().genId())),
      _momentum_threshold(config().momentum_threshold()){
    /**/
  }

  // return true if particle is associated with the configured generator
  bool GenIdSimParticleSelectionTool::Select(const SimParticle& particle){
    bool rv = false;
    const auto& gp = particle.genParticle();
    if (gp.isNonnull()){
      rv = (gp->generatorId() == _genId);
    }
    return rv;
  }
} // namespace mu2e

DEFINE_ART_CLASS_TOOL(mu2e::GenIdSimParticleSelectionTool)
