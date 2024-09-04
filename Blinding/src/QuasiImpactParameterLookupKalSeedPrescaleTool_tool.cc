// Ed Callaghan
// art tool to prescale KalSeeds via a lookup-table on the ~minimum DS radial coordinate
// September 2024

#include "Offline/Blinding/inc/QuasiImpactParameterLookupKalSeedPrescaleTool.hh"

namespace mu2e{
  QuasiImpactParameterLookupKalSeedPrescaleTool::QuasiImpactParameterLookupKalSeedPrescaleTool(const Parameters& config):
      UnivariateLookupKalSeedPrescaleTool(config().radius_start(),
                                          config().radius_step_size(),
                                          config().prescale()){
    for (const auto& name: config().surface_ids()){
      _surface_ids.emplace_back(name);
    }
  }

  double QuasiImpactParameterLookupKalSeedPrescaleTool::calculate_observable(const KalSeed& kalseed){
    KinKal::LoopHelix helix;
    const auto& end = kalseed.intersections().end();
    auto sit = _surface_ids.begin();
    bool adequate = false;
    while ((!adequate) && (sit != _surface_ids.end())){
      auto const& iit = kalseed.intersection(*sit);
      if (iit != end){
        helix = iit->loopHelix();
        adequate = true;
      }
      sit++;
    }

    // if no intersectiojn found, then force-zero the upstream rv
    auto rv = 0.1 * _xmin;
    if (adequate){
      rv = hypot(helix.cx(), helix.cy()) - fabs(helix.rad());
    }
    return rv;
  }
} // namespace mu2e

DEFINE_ART_CLASS_TOOL(mu2e::QuasiImpactParameterLookupKalSeedPrescaleTool)
