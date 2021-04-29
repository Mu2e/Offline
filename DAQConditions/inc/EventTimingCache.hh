#ifndef DAQConditions_EventTimingCache_hh
#define DAQConditions_EventTimingCache_hh

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
#include "DAQConditions/inc/EventTimingMaker.hh"


namespace mu2e {
  class EventTimingCache : public ProditionsCache {
  public: 
    EventTimingCache(EventTimingConfig const& config):
      ProditionsCache(EventTiming::cxname,config.verbose()),
      _useDb(config.useDb()),_maker(config) {}

    void initialize() {
    }

    set_t makeSet(art::EventID const& eid) {
      return ProditionsEntity::set_t();
    }

    DbIoV makeIov(art::EventID const& eid) {
      DbIoV iov;
      iov.setMax();
      return iov;
    }

    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      return _maker.fromFcl();
    }

  private:
    bool _useDb;
    EventTimingMaker _maker;
  };
};

#endif
