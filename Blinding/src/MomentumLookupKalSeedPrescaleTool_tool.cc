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
    // search for downward-going KalIntersection with a matching surface
    KalIntersection intersection;
    auto sit = _surface_ids.begin();
    bool adequate = false;
    while ((!adequate) && (sit != _surface_ids.end())){
      const auto& iits = kalseed.intersections(*sit);
      // there are at most 2 valid intersections (upgoing / downgoing)
      // so selecting on which one is downgoing is sufficient
      for (const auto& iit: iits){
        auto direction = iit->momentum3().Unit().Z();
        if (0 < direction){
          intersection = *iit;
          adequate = true;
        }
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
