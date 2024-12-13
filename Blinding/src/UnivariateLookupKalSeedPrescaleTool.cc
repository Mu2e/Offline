// Ed Callaghan
// Interface for art tool to calculate a single-variable functional prescale
// September 2024

#include "Offline/Blinding/inc/UnivariateLookupKalSeedPrescaleTool.hh"

namespace mu2e{
  // TODO
  // read in from GlobalConstantsService, or third-party
  // propagating from config is just ridiculous
  UnivariateLookupKalSeedPrescaleTool::UnivariateLookupKalSeedPrescaleTool(double xmin, double step_size, std::vector<double> coordinates){
    _xmin = xmin;
    _xmax = _xmin + step_size * coordinates.size();
    // can throw
    _spline = cardinal_cubic_b_spline(coordinates.begin(), coordinates.end(),
                                      xmin, step_size);
  }

  double UnivariateLookupKalSeedPrescaleTool::calculate_acceptance_rate(const KalSeed& kalseed){
    double x = this->calculate_observable(kalseed);
    double rv = 0.0;
    if ((_xmin <= x) && (x < _xmax)){
      rv = _spline(x);
    }
    return rv;
  }
} // namespace mu2e
