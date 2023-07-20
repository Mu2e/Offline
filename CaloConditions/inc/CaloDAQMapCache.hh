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
    }

    set_t makeSet(art::EventID const& eid) {
      ProditionsEntity::set_t cids;
      return cids;
    }

    DbIoV makeIov(art::EventID const& eid) {
      DbIoV iov;
      iov.setMax(); // start with full IOV range
      return iov;
    }

    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      if (_useDb) {
        return _maker.fromDb();
      } else {
        return _maker.fromFcl();
      }
    }

  private:
    bool _useDb;
    CaloDAQMapMaker _maker;

  };
}

#endif
