#ifndef TrackerConditions_TrkPanelMapMaker_hh
#define TrackerConditions_TrkPanelMapMaker_hh

//
// construct a TrkPanelMap conditions entity
// from fcl or database
//

#include "Offline/DbTables/inc/TrkPanelMap.hh"
#include "Offline/TrackerConfig/inc/TrkPanelMapConfig.hh"
#include "Offline/TrackerConditions/inc/TrkPanelMapEntity.hh"

namespace mu2e {

  class TrkPanelMapMaker {
  public:
    TrkPanelMapMaker(TrkPanelMapConfig const& config):config_(config) {}
    TrkPanelMapMaker() {}
    TrkPanelMapEntity::ptr_t fromFcl();
    TrkPanelMapEntity::ptr_t fromDb (TrkPanelMap::cptr_t Table);

  private:
      const TrkPanelMapConfig config_;
    
  };
}

#endif

