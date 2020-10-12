#ifndef SimulationConditions_SimBookkeeperCache_hh
#define SimulationConditions_SimBookkeeperCache_hh

//
// SimBookkeeperCache for ProditionsCache
//
#include <iostream>

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
#include "DbTables/inc/DbIoV.hh"
#include "DbService/inc/DbHandle.hh"
#include "DbTables/inc/SimEfficiencies.hh"
#include "SimulationConditions/inc/SimBookkeeperMaker.hh"

namespace mu2e {

  class SimBookkeeperCache : public ProditionsCache {
  public:
    SimBookkeeperCache(SimBookkeeperConfig const& config):
      ProditionsCache("SimBookkeeper", config.verbose()),
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
    std::unique_ptr<DbHandle<SimEfficiencies> > _tqDb_p;
  };
};

#endif
