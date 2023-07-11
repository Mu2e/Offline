#ifndef TrackerConditions_StrawElectronicsCache_hh
#define TrackerConditions_StrawElectronicsCache_hh

#include "Offline/Mu2eInterfaces/inc/ProditionsCache.hh"
#include "Offline/DbTables/inc/DbIoV.hh"
#include "Offline/DbService/inc/DbHandle.hh"
#include "Offline/DbTables/inc/TrkDelayPanel.hh"
#include "Offline/DbTables/inc/TrkDelayRStraw.hh"
#include "Offline/DbTables/inc/TrkPreampStraw.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"

#include "Offline/TrackerConditions/inc/StrawElectronicsMaker.hh"


namespace mu2e {
  class StrawElectronicsCache : public ProditionsCache {
    public:
      StrawElectronicsCache(StrawElectronicsConfig const& config):
        ProditionsCache(StrawElectronics::cxname,config.verbose()),
        _useDb(config.useDb()),_maker(config) {}

      void initialize() {
        _eventTiming_p = std::make_unique<ProditionsHandle<EventTiming> >();
        if(_useDb) {
          _tdp_p  = std::make_unique<DbHandle<TrkDelayPanel>>();
          _tdrs_p = std::make_unique<DbHandle<TrkDelayRStraw>>();
          _tps_p  = std::make_unique<DbHandle<TrkPreampStraw>>();
        }
      }

      set_t makeSet(art::EventID const& eid) {
        auto et = _eventTiming_p->get(eid);
        ProditionsEntity::set_t cids = et.getCids();
        if(_useDb) { // use fcl config, overwrite part from DB
          // get the tables up to date
          _tdp_p->get(eid);
          _tdrs_p->get(eid);
          _tps_p->get(eid);
          // save which data goes into this instance of the service
          cids.insert(_tdp_p->cid());
          cids.insert(_tdrs_p->cid());
          cids.insert(_tps_p->cid());
        }
        return cids;
      }

      DbIoV makeIov(art::EventID const& eid) {
        _eventTiming_p->get(eid);
        DbIoV iov = _eventTiming_p->iov();
        if(_useDb) { // use fcl config, overwrite part from DB
          // get the tables up to date
          _tdp_p->get(eid);
          _tdrs_p->get(eid);
          _tps_p->get(eid);
          // restrict the valid range ot the overlap
          iov.overlap(_tdp_p->iov());
          iov.overlap(_tdrs_p->iov());
          iov.overlap(_tps_p->iov());
        }
        return iov;
      }

      ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
        auto et = _eventTiming_p->getPtr(eid);
        if(_useDb) {
          return _maker.fromDb( _tdp_p->getPtr(eid),
              _tdrs_p->getPtr(eid),
              _tps_p->getPtr(eid),
              et);
        } else {
          return _maker.fromFcl(et);
        }
      }

    private:
      bool _useDb;
      StrawElectronicsMaker _maker;


      // these handles are not default constructed
      // so to not create a dependency loop on construction
      std::unique_ptr<ProditionsHandle<EventTiming> > _eventTiming_p;
      // these handles are not default constructed
      // so the db can be completely turned off
      std::unique_ptr<DbHandle<TrkDelayPanel>> _tdp_p;
      std::unique_ptr<DbHandle<TrkDelayRStraw>> _tdrs_p;
      std::unique_ptr<DbHandle<TrkPreampStraw>> _tps_p;

  };
}

#endif
