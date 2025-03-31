// Ed Callaghan
// Produce sinusoidal waves
// February 2025

#include "Offline/TrackerMC/inc/TruncatedSinusoidTool.hh"

namespace mu2e{
  TruncatedSinusoidTool::TruncatedSinusoidTool(const Parameters& config):
      _amplitude(config().amplitude()),
      _frequency(config().frequency()),
      _phase(config().phase()),
      _time_lo(config().time_lo()),
      _time_hi(config().time_hi()){
    /**/
  }

  UnaryFunctionPtr TruncatedSinusoidTool::Sample(){
    auto rv = std::make_shared<TruncatedSinusoid>(_amplitude,_frequency, _phase,
                                                  _time_lo, _time_hi);
    return rv;
  }
} // namespace mu2e

DEFINE_ART_CLASS_TOOL(mu2e::TruncatedSinusoidTool)
