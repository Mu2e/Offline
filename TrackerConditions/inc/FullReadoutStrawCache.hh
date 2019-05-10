#ifndef TrackerConditions_FullReadoutStrawCache_hh
#define TrackerConditions_FullReadoutStrawCache_hh

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
//#include "DbService/inc/DbHandle.hh"
#include "TrackerConditions/inc/FullReadoutStrawMaker.hh"


namespace mu2e {
  class FullReadoutStrawCache : public ProditionsCache {
  public: 
    FullReadoutStrawCache(FullReadoutStrawConfig const& config):
      ProditionsCache("FullReadoutStraw",config.verbose()),
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
    FullReadoutStrawMaker _maker;

  };
};

#endif
