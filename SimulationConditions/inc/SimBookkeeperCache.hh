#ifndef SimulationConditions_SimBookkeeperCache_hh
#define SimulationConditions_SimBookkeeperCache_hh

//
// SimBookkeeperCache for ProditionsCache
//

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

    void initialize() {
      if(_useDb) {
        _tqDb_p  = std::make_unique<DbHandle<SimEfficiencies> >();
      }
    }

    set_t makeSet(art::EventID const& eid) {
      ProditionsEntity::set_t cids;
      if(_useDb) { // use fcl config, overwrite part from DB
        // get the tables up to date
        _tqDb_p->get(eid);
        // save which data goes into this instance of the service
        cids.insert(_tqDb_p->cid());
      }
      return cids;
    }

    DbIoV makeIov(art::EventID const& eid) {
      DbIoV iov;
      iov.setMax(); // start with full IOV range
      if(_useDb) { // use fcl config, overwrite part from DB
        // get the tables up to date
        _tqDb_p->get(eid);
        // restrict the valid range ot the overlap
        iov.overlap(_tqDb_p->iov());
      }
      return iov;
    }

    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      if(_useDb) {
        return _maker.fromDb( _tqDb_p->getPtr(eid) );
      } else {
        return _maker.fromFcl();
      }
    }

  private:
    bool _useDb;
    SimBookkeeperMaker _maker;

    // these handles are not default constructed
    // so the db can be completely turned off
    std::unique_ptr<DbHandle<SimEfficiencies> > _tqDb_p;
  };
};

#endif
