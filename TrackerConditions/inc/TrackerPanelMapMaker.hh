#ifndef TrackerConditions_TrackerPanelMapMaker_hh
#define TrackerConditions_TrackerPanelMapMaker_hh

//
// construct a TrackerPanelMap conditions entity
// from fcl or database
//

#include "Offline/DbTables/inc/TrkPanelMap.hh"
#include "Offline/TrackerConfig/inc/TrackerPanelMapConfig.hh"
#include "Offline/TrackerConditions/inc/TrackerPanelMap.hh"

namespace mu2e {

  class TrackerPanelMapMaker {
  public:
    TrackerPanelMapMaker(TrackerPanelMapConfig const& config):config_(config) {}
    TrackerPanelMapMaker() {}
    TrackerPanelMap::ptr_t fromFcl();
    TrackerPanelMap::ptr_t fromDb (TrkPanelMap::cptr_t Table);

  private:
      const TrackerPanelMapConfig config_;

  };
}

#endif

