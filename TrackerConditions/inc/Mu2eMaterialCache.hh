#ifndef TrackerConditions_Mu2eMaterialCache_hh
#define TrackerConditions_Mu2eMaterialCache_hh

//
// This Proditions entitiy cache should only ever hold one
// because Mu2eMaterial has pointers to BTrk singletons
//

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
#include "TrackerConditions/inc/Mu2eMaterialMaker.hh"


namespace mu2e {
  class Mu2eMaterialCache : public ProditionsCache {
  public: 
    Mu2eMaterialCache(Mu2eMaterialConfig const& config):
      ProditionsCache("Mu2eMaterial",config.verbose()),
      _maker(config) {
      // force a fake update so BTrk TrkParticle 
      // is working when modules are created
      update(art::EventID(1,0,0));
    }
    
    void initialize() {
    }
    set_t makeSet(art::EventID const& eid) {
      return ProditionsEntity::set_t();
    }
    DbIoV makeIov(art::EventID const& eid) {
      DbIoV iov;
      iov.setMax(); // all runs
      return iov;
    }
    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      return _maker.fromFcl();
    }

  private:
    Mu2eMaterialMaker _maker;

  };
};

#endif
