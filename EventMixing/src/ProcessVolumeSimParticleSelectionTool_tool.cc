// Ed Callaghan
// Tool to select SimParticles created by a specific G4 physics process
// November 2024

#include "Offline/EventMixing/inc/ProcessVolumeSimParticleSelectionTool.hh"

namespace mu2e{

  ProcessVolumeSimParticleSelectionTool::ProcessVolumeSimParticleSelectionTool(const Parameters& config):
      _processCode(ProcessCode::findByName(config().process())),
      _volume(config().volume()),
      _momentum_threshold(config().momentum_threshold()){
    // internally-managed hard config for temporary functionality
    fhicl::ParameterSet subconfig;
    subconfig.put("tool_type", "InnerProtonAbsorberPseudoVolumeLookupTool");
    subconfig.put("IPA", "protonabs1");
    // beginning of sha512(protonabs1), to not clobber any reasonable use-case
    subconfig.put("Other", "e536f774a6");
    _lookup = art::make_tool<InnerProtonAbsorberPseudoVolumeLookupTool>(subconfig);
  }

  // return true if particle matches the configured start process and volume
  // the queries are nested to condition the more-expensive coordinate
  // manipulations
  bool ProcessVolumeSimParticleSelectionTool::Select(const SimParticle& particle){
    bool rv = false;
    bool process_matched = (particle.creationCode() == this->_processCode);
    if (process_matched){
      double momentum = particle.startMomentum().vect().mag();
      bool momentum_passed = (_momentum_threshold < momentum);
      if (momentum_passed){
        bool volume_matched = (_lookup->StartVolume(particle) == this->_volume);
        if (volume_matched){
          rv = true;
        }
      }
    }
    return rv;
  }
} // namespace mu2e

DEFINE_ART_CLASS_TOOL(mu2e::ProcessVolumeSimParticleSelectionTool)
