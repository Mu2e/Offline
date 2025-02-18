// Ed Callaghan
// Produce sinusoidal waves
// February 2025

#include "Offline/TrackerMC/inc/SinusoidalWireSignalTool.hh"

namespace mu2e{
  SinusoidalWireSignalTool::SinusoidalWireSignalTool(const Parameters& config):
      _amplitude(config().amplitude()),
      _frequency(config().frequency()),
      _time_lo(config().time_lo()),
      _time_hi(config().time_hi()){
    /**/
  }

  AnalogWireSignalPtr SinusoidalWireSignalTool::Sample(){
    AnalogWireSignalPtr rv = std::make_shared<SinusoidalWireSignal>(_amplitude, _frequency, 0.0, _time_lo, _time_hi);
    return rv;
  }
} // namespace mu2e

DEFINE_ART_CLASS_TOOL(mu2e::SinusoidalWireSignalTool)
