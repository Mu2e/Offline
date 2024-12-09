// Ed Callaghan
// Select all detector steps
// November 2024

#include "Offline/EventMixing/inc/UniversalDetectorStepSelectionTool.hh"

namespace mu2e{
  UniversalDetectorStepSelectionTool::UniversalDetectorStepSelectionTool(const Parameters& config){
    /**/
  }

  bool UniversalDetectorStepSelectionTool::Select(const CaloShowerStep& step){
    bool rv = this->select(step);
    return rv;
  }

  bool UniversalDetectorStepSelectionTool::Select(const CrvStep& step){
    bool rv = this->select(step);
    return rv;
  }

  bool UniversalDetectorStepSelectionTool::Select(const StrawGasStep& step){
    bool rv = this->select(step);
    return rv;
  }
}

DEFINE_ART_CLASS_TOOL(mu2e::UniversalDetectorStepSelectionTool)
