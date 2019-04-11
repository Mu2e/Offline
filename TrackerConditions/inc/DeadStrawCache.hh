#ifndef TrackerConditions_DeadStrawCache_hh
#define TrackerConditions_DeadStrawCache_hh

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
//#include "DbService/inc/DbHandle.hh"
#include "TrackerConditions/inc/DeadStrawMaker.hh"


namespace mu2e {
  class DeadStrawCache : public ProditionsCache {
  public: 
    DeadStrawCache(DeadStrawConfig const& config):
      ProditionsCache("DeadStraw",config.verbose()),
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
    DeadStrawMaker _maker;

  };
};

#endif
