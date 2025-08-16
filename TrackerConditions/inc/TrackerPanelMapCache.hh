#ifndef TrackerConditions_TrackerPanelMapCache_hh
#define TrackerConditions_TrackerPanelMapCache_hh

#include "Offline/Mu2eInterfaces/inc/ProditionsCache.hh"
#include "Offline/DbTables/inc/DbIoV.hh"
#include "Offline/DbService/inc/DbHandle.hh"

#include "Offline/TrackerConditions/inc/TrackerPanelMapMaker.hh"
#include "Offline/TrackerConfig/inc/TrackerPanelMapConfig.hh"


namespace mu2e {
  class TrackerPanelMapCache : public ProditionsCache {
    public:
      TrackerPanelMapCache(TrackerPanelMapConfig const& config):
        ProditionsCache(TrackerPanelMap::cxname,config.verbose()),
        _useDb(config.useDb()),_maker(config) {}

      void initialize() {
        if(_useDb) {
          _tpm_p  = std::make_unique<DbHandle<TrkPanelMap>>();
        }
      }

      set_t makeSet(art::EventID const& eid) {
        ProditionsEntity::set_t cids;
        if(_useDb) {                    // use fcl config, overwrite part from DB
                                        // get the tables up to date
          _tpm_p->get(eid);
                                        // save which data goes into this instance of the service
          cids.insert(_tpm_p->cid());
        }
        return cids;
      }

      DbIoV makeIov(art::EventID const& eid) {
        DbIoV iov;
        if(_useDb) { // use fcl config, overwrite part from DB
          // get the tables up to date
          _tpm_p->get(eid);
          iov.overlap(_tpm_p->iov());
        }
        return iov;
      }

      ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
        if(_useDb) {
          return _maker.fromDb( _tpm_p->getPtr(eid));
        }
        else {
          return _maker.fromFcl();
        }
      }

    private:
      bool                 _useDb;
      TrackerPanelMapMaker _maker;

      // these handles are not default constructed
      // so the db can be completely turned off
      std::unique_ptr<DbHandle<TrkPanelMap>> _tpm_p;
  };
}

#endif
