// Ed Callaghan
// art tool to prescale KalSeeds via a lookup-table on the momentum
// September 2024

#include "Offline/Blinding/inc/MomentumLookupKalSeedPrescaleTool.hh"

namespace mu2e{
  MomentumLookupKalSeedPrescaleTool::MomentumLookupKalSeedPrescaleTool(const Parameters& config):
      UnivariateLookupKalSeedPrescaleTool(config().momentum_start(),
                                          config().momentum_step_size(),
                                          config().prescale()){
    for (const auto& name: config().surface_ids()){
      _surface_ids.emplace_back(name);
    }
  }

  double MomentumLookupKalSeedPrescaleTool::calculate_observable(const KalSeed& kalseed){
    KalIntersection intersection;
    const auto& end = kalseed.intersections().end();
    auto sit = _surface_ids.begin();
    bool adequate = false;
    while ((!adequate) && (sit != _surface_ids.end())){
      auto const& iit = kalseed.intersection(*sit);
      if (iit != end){
        intersection = *iit;
        adequate = true;
      }
      sit++;
    }

    // if no intersectiojn found, then force-zero the upstream rv
    auto rv = 0.1 * _xmin;
    if (adequate){
      rv = intersection.mom();
    }
    return rv;
  }
} // namespace mu2e

DEFINE_ART_CLASS_TOOL(mu2e::MomentumLookupKalSeedPrescaleTool)
