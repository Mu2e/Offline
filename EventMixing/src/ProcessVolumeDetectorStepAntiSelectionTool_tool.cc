// Ed Callaghan
// Select detector steps which are not downstream of a certain process
// November 2024

#include "Offline/EventMixing/inc/ProcessVolumeDetectorStepAntiSelectionTool.hh"

namespace mu2e{
  ProcessVolumeDetectorStepAntiSelectionTool::ProcessVolumeDetectorStepAntiSelectionTool(const Parameters& config):
      _energy_lo(config().energy_lo()),
      _energy_hi(config().energy_hi()){
    for (const auto& processCode: config().processes()){
      _processCodes.insert(ProcessCode::findByName(processCode).id());
    }
    for (const auto& volume: config().volumes()){
      _volumes.insert(volume);
    }

    // internally-managed hard config for temporary functionality
    fhicl::ParameterSet subconfig;
    subconfig.put("tool_type", "PseudoCylindricalVolumeLookupTool");
    subconfig.put("IPA", "protonabs1");
    subconfig.put("ST", "Foil_00");
    // beginning of sha512(protonabs1), to not clobber any reasonable use-case
    subconfig.put("Other", "e536f774a6");
    _lookup = art::make_tool<PseudoCylindricalVolumeLookupTool>(subconfig);
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
    const auto& processCode = particle.creationCode();
    bool process_matched = (0 < _processCodes.count(processCode.id()));
    if (process_matched){
      double energy = particle.startMomentum().e();
      bool energy_passed = true;
      if (_energy_lo.has_value()){
        energy_passed = energy_passed && (_energy_lo.value() < energy);
      }
      if (_energy_hi.has_value()){
        energy_passed = energy_passed && (energy < _energy_hi.value());
      }
      if (energy_passed){
        const auto volume = _lookup->StartVolume(particle);
        bool volume_matched = (0 < _volumes.count(volume));
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
