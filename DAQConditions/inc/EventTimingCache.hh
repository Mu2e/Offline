#ifndef DAQConditions_EventTimingCache_hh
#define DAQConditions_EventTimingCache_hh

#include "Offline/Mu2eInterfaces/inc/ProditionsCache.hh"
#include "Offline/DAQConditions/inc/EventTimingMaker.hh"
#include "Offline/DbTables/inc/DAQTiming.hh"


namespace mu2e {
  class EventTimingCache : public ProditionsCache {
  public:
    EventTimingCache(EventTimingConfig const& config):
      ProditionsCache(EventTiming::cxname,config.verbose()),
      _useDb(config.useDb()),_maker(config) {}

    void initialize() {
      if (_useDb) {
        _dt_p = std::make_unique<DbHandle<DAQTiming>>();
      }
    }

    set_t makeSet(art::EventID const& eid) {
      ProditionsEntity::set_t cids;
      if (_useDb) {
        _dt_p->get(eid);
        cids.insert(_dt_p->cid());
      }
      return cids;
    }

    DbIoV makeIov(art::EventID const& eid) {
      DbIoV iov;
      iov.setMax();
      if (_useDb) {
        _dt_p->get(eid);
        iov.overlap(_dt_p->iov());
      }
      return iov;
    }

    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      if (_useDb) {
        return _maker.fromDb(_dt_p->getPtr(eid));
      } else {
        return _maker.fromFcl();
      }
    }

  private:
    bool _useDb;
    EventTimingMaker _maker;
    std::unique_ptr<DbHandle<DAQTiming>> _dt_p;
  };
}

#endif
