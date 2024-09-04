// Ed Callaghan
// art tool to prescale KalSeeds via a lookup-table on the tracker exit azimuthal coordinate
// September 2024

#include "Offline/Blinding/inc/TrackerExitAzimuthLookupKalSeedPrescaleTool.hh"

namespace mu2e{
  TrackerExitAzimuthLookupKalSeedPrescaleTool::TrackerExitAzimuthLookupKalSeedPrescaleTool(const Parameters& config):
      UnivariateLookupKalSeedPrescaleTool(config().azimuth_start(),
                                          config().azimuth_step_size(),
                                          config().prescale()){
    for (const auto& name: config().surface_ids()){
      _surface_ids.emplace_back(name);
    }
  }

  double TrackerExitAzimuthLookupKalSeedPrescaleTool::calculate_observable(const KalSeed& kalseed){
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
      // TODO explicitly cast to detector coordinates, instead of the
      // hacked-in rotation used here
      auto position = intersection.position3().Unit();
      rv = atan2(position.x(), position.y());
    }
    return rv;
  }
} // namespace mu2e

DEFINE_ART_CLASS_TOOL(mu2e::TrackerExitAzimuthLookupKalSeedPrescaleTool)
