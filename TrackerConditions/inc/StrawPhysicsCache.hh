#ifndef TrackerConditions_StrawPhysicsCache_hh
#define TrackerConditions_StrawPhysicsCache_hh

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"
#include "TrackerConditions/inc/StrawDrift.hh"
#include "TrackerConditions/inc/StrawPhysicsMaker.hh"


namespace mu2e {
  class StrawPhysicsCache : public ProditionsCache {
  public: 
    StrawPhysicsCache(StrawPhysicsConfig const& config):
      ProditionsCache("StrawPhysics",config.verbose()),
      _useDb(config.useDb()),_maker(config) {}

    //    std::string const& name() const { return _name; }


    void initialize() {
      _strawDrift_p = std::make_unique<ProditionsHandle<StrawDrift> >();
    }
    set_t makeSet(art::EventID const& eid) {
     return _strawDrift_p->get(eid).getCids();
    }
    DbIoV makeIov(art::EventID const& eid) {
      _strawDrift_p->get(eid); // make sure it is current
      return _strawDrift_p->iov();
    }
    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      auto sd = _strawDrift_p->getPtr(eid);
      return _maker.fromFcl(sd);
    }

  private:
    bool _useDb;
    StrawPhysicsMaker _maker;

    // this handle is not default constructed
    // so to not create a dependency loop on construction
    std::unique_ptr<ProditionsHandle<StrawDrift> > _strawDrift_p;
  };
};

#endif
