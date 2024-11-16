// Ed Callaghan
// Select detector steps which are not downstream of a certain process
// November 2024

#include "Offline/EventMixing/inc/DescendantDetectorStepAntiSelectionTool.hh"

namespace mu2e{
  DescendantDetectorStepAntiSelectionTool::DescendantDetectorStepAntiSelectionTool(const Parameters& config){
    auto subconfig = config().selection.get<fhicl::ParameterSet>();
    _selection = art::make_tool<SimParticleSelectionTool>(subconfig);
  }

  bool DescendantDetectorStepAntiSelectionTool::Select(const CaloShowerStep& step){
    bool rv = this->select(step);
    return rv;
  }

  bool DescendantDetectorStepAntiSelectionTool::Select(const CrvStep& step){
    bool rv = this->select(step);
    return rv;
  }

  bool DescendantDetectorStepAntiSelectionTool::Select(const StrawGasStep& step){
    bool rv = this->select(step);
    return rv;
  }

  bool DescendantDetectorStepAntiSelectionTool::particle_match(const SimParticle& particle){
    bool rv = _selection->Select(particle);
    return rv;
  }

  bool DescendantDetectorStepAntiSelectionTool::heritage_match(const SimParticle& particle){
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

DEFINE_ART_CLASS_TOOL(mu2e::DescendantDetectorStepAntiSelectionTool)
