#ifndef CaloConditions_CaloDAQMapCache_hh
#define CaloConditions_CaloDAQMapCache_hh

#include "Offline/Mu2eInterfaces/inc/ProditionsCache.hh"
//#include "DbTables/inc/DbIoV.hh"
#include "Offline/DbService/inc/DbHandle.hh"
#include "Offline/DbTables/inc/CalRoIDMapDIRACToOffline.hh"
#include "Offline/DbTables/inc/CalRoIDMapOfflineToDIRAC.hh"

#include "Offline/CaloConditions/inc/CaloDAQMapMaker.hh"


namespace mu2e {
  class CaloDAQMapCache : public ProditionsCache {
  public:
    CaloDAQMapCache(CaloDAQMapConfig const& config):
    ProditionsCache(CaloDAQMap::cxname,config.verbose()),
      _useDb(config.useDb()),_maker(config) {}

    void initialize() {
      if(_useDb) {
        _tdtc_p = std::make_unique<DbHandle<CalRoIDMapDIRACToOffline>>();
        _tctd_p = std::make_unique<DbHandle<CalRoIDMapOfflineToDIRAC>>();
      }
    }

    set_t makeSet(art::EventID const& eid) {
      ProditionsEntity::set_t cids;
      if(_useDb) { // use fcl config, overwrite part from DB
        // get the tables up to date
        _tdtc_p->get(eid);
        _tctd_p->get(eid);
        // save which data goes into this instance of the service
        cids.insert(_tdtc_p->cid());
        cids.insert(_tctd_p->cid());
      }
      return cids;
    }

    DbIoV makeIov(art::EventID const& eid) {
      DbIoV iov;
      iov.setMax(); // start with full IOV range
      if(_useDb) { // use fcl config, overwrite part from DB
        // get the tables up to date
        _tdtc_p->get(eid);
        _tctd_p->get(eid);
        // restrict the valid range ot the overlap
        iov.overlap(_tdtc_p->iov());
        iov.overlap(_tctd_p->iov());
      }
      return iov;
    }

    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      if(_useDb) {
        return _maker.fromDb( _tdtc_p->getPtr(eid),
                              _tctd_p->getPtr(eid) );
      } else {
        return _maker.fromFcl();
      }
    }

  private:
    bool _useDb;
    CaloDAQMapMaker _maker;

    // these handles are not default constructed
    // so the db can be completely turned off
    std::unique_ptr<DbHandle<CalRoIDMapDIRACToOffline>> _tdtc_p;
    std::unique_ptr<DbHandle<CalRoIDMapOfflineToDIRAC>> _tctd_p;

  };
}

#endif
