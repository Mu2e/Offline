// Ed Callaghan
// Select detector steps which are not downstream of a certain process
// November 2024

#include "Offline/EventMixing/inc/ProcessVolumeDetectorStepAntiSelectionTool.hh"

namespace mu2e{
  ProcessVolumeDetectorStepAntiSelectionTool::ProcessVolumeDetectorStepAntiSelectionTool(const Parameters& config):
      _processCode(ProcessCode::findByName(config().process())),
      _volume(config().volume()),
      _momentum_threshold(config().momentum_threshold()){
    //
    fhicl::ParameterSet subconfig;
    subconfig.put("tool_type", "InnerProtonAbsorberPseudoVolumeLookupTool");
    subconfig.put("IPA", "protonabs1");
    // beginning of sha512(protonabs1), to not clobber any reasonable use-case
    subconfig.put("Other", "e536f774a6");
    _lookup = art::make_tool<InnerProtonAbsorberPseudoVolumeLookupTool>(subconfig);
  }

  bool ProcessVolumeDetectorStepAntiSelectionTool::Select(const CaloShowerStep& step){
    bool rv = this->select(step);
    return rv;
  }

  bool ProcessVolumeDetectorStepAntiSelectionTool::Select(const CrvStep& step){
    bool rv = this->select(step);
    return rv;
  }

  bool ProcessVolumeDetectorStepAntiSelectionTool::Select(const StrawGasStep& step){
    bool rv = this->select(step);
    return rv;
  }

  // return true if particle matches the configured start process and volume
  // the queries are nested to condition the more-expensive coordinate
  // manipulations
  bool ProcessVolumeDetectorStepAntiSelectionTool::particle_match(const SimParticle& particle){
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

  bool ProcessVolumeDetectorStepAntiSelectionTool::heritage_match(const SimParticle& particle){
    // check condition of current particle
    if (this->particle_match(particle)){
      return true;
    }

    // if not a match, defer match to its parent, if available
    const auto& parent = particle.parent();
    if (parent.isNull()){
      return false;
    }
    else{
      return this->heritage_match(*parent);
    }
  }
} // namespace

DEFINE_ART_CLASS_TOOL(mu2e::ProcessVolumeDetectorStepAntiSelectionTool)
