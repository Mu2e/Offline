//
// SimBookkeeperCache for ProditionsCache
//

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "SimulationConditions/inc/SimBookkeeperCache.hh"

namespace mu2e {

  void SimBookkeeperCache::initialize() {
    if(_useDb) {
      _tqDb_p  = std::make_unique<DbHandle<SimEfficiencies> >();
    }
  }

  ProditionsEntity::set_t SimBookkeeperCache::makeSet(art::EventID const& eid) {
    ProditionsEntity::set_t cids;
    if(_useDb) { // use fcl config, overwrite part from DB
      // get the tables up to date
      _tqDb_p->get(eid);
      // save which data goes into this instance of the service
      cids.insert(_tqDb_p->cid());
    }
    return cids;
  }

  DbIoV SimBookkeeperCache::makeIov(art::EventID const& eid) {
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

  ProditionsEntity::ptr SimBookkeeperCache::makeEntity(art::EventID const& eid) {
    if(_useDb) {
      return _maker.fromDb( _tqDb_p->getPtr(eid) );
    } else {
      return _maker.fromFcl();
    }
  }
}
