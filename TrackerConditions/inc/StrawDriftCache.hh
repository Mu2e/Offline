#ifndef TrackerConditions_StrawDriftCache_hh
#define TrackerConditions_StrawDriftCache_hh

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
//#include "DbService/inc/DbHandle.hh"
#include "TrackerConditions/inc/StrawDriftMaker.hh"


namespace mu2e {
  class StrawDriftCache : public ProditionsCache {
  public: 
    StrawDriftCache(StrawDriftConfig const& config):
      ProditionsCache("StrawDrift",config.verbose()),
      _useDb(config.useDb()),_maker(config) {}


    void initialize() {
    }
    set_t makeSet(art::EventID const& eid) {
      return set_t();
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
    StrawDriftMaker _maker;

  };
};

#endif
