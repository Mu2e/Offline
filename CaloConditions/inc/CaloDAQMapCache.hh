#ifndef CaloConditions_CaloDAQMapCache_hh
#define CaloConditions_CaloDAQMapCache_hh

#include "Offline/Mu2eInterfaces/inc/ProditionsCache.hh"
//#include "DbTables/inc/DbIoV.hh"
#include "Offline/DbService/inc/DbHandle.hh"

#include "Offline/CaloConditions/inc/CaloDAQMapMaker.hh"


namespace mu2e {
  class CaloDAQMapCache : public ProditionsCache {
  public:
    CaloDAQMapCache(CaloDAQMapConfig const& config):
    ProditionsCache(CaloDAQMap::cxname,config.verbose()),
      _useDb(config.useDb()),_maker(config) {}

    void initialize() {
      if (_useDb) {
        _cch_p = std::make_unique<DbHandle<CalChannels>>();
      }
    }

    set_t makeSet(art::EventID const& eid) {
      ProditionsEntity::set_t cids;
      if (_useDb) {
        _cch_p->get(eid);
        cids.insert(_cch_p->cid());
      }
     return cids;
    }

    DbIoV makeIov(art::EventID const& eid) {
      DbIoV iov;
      iov.setMax(); // start with full IOV range
      if (_useDb) {
        _cch_p->get(eid);
        iov.overlap(_cch_p->iov());
      }
      return iov;
    }

    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      if (_useDb) {
        return _maker.fromDb(_cch_p->getPtr(eid));
      } else {
        return _maker.fromFcl();
      }
    }

  private:
    bool _useDb;
    CaloDAQMapMaker _maker;

    // these handles are not default constructed
    // so the db can be completely turned off
    std::unique_ptr<DbHandle<CalChannels>> _cch_p;

  };
}

#endif
