#ifndef SimulationConditions_SimBookkeeperCache_hh
#define SimulationConditions_SimBookkeeperCache_hh

//
// Caches the SimBookeeper ProditionsEntity
//
#include <iostream>

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "Offline/Mu2eInterfaces/inc/ProditionsCache.hh"
#include "Offline/DbTables/inc/DbIoV.hh"
#include "Offline/DbService/inc/DbHandle.hh"
#include "Offline/DbTables/inc/SimEfficiencies2.hh"
#include "Offline/SimulationConditions/inc/SimBookkeeperMaker.hh"

namespace mu2e {

  class SimBookkeeperCache : public ProditionsCache {
  public:
    SimBookkeeperCache(SimBookkeeperConfig const& config):
      ProditionsCache(SimBookkeeper::cxname, config.verbose()),
      _useDb(config.useDb()),_maker(config) {}

    void initialize();
    ProditionsEntity::set_t makeSet(art::EventID const& eid);
    DbIoV makeIov(art::EventID const& eid);
    ProditionsEntity::ptr makeEntity(art::EventID const& eid);

  private:
    bool _useDb;
    SimBookkeeperMaker _maker;

    // these handles are not default constructed
    // so the db can be completely turned off
    std::unique_ptr<DbHandle<SimEfficiencies2> > _tqDb_p;
  };
}

#endif
