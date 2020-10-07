#ifndef TrackerConditions_TrackerStatusCache_hh
#define TrackerConditions_TrackerStatusCache_hh

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
#include "TrackerConditions/inc/TrackerStatusMaker.hh"


namespace mu2e {
  class TrackerStatusCache : public ProditionsCache {
  public: 
    TrackerStatusCache(TrackerStatusConfig const& config):
      ProditionsCache("TrackerStatus",config.settings().verbose()),
      _useDb(config.settings().useDb()),_maker(config) {}


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
    TrackerStatusMaker _maker;

  };
};

#endif
