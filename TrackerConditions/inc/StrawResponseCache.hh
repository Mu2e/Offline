#ifndef TrackerConditions_StrawResponseCache_hh
#define TrackerConditions_StrawResponseCache_hh

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
#include "TrackerConditions/inc/StrawResponseMaker.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"


namespace mu2e {
  class StrawResponseCache : public ProditionsCache {
  public: 
    StrawResponseCache(StrawResponseConfig const& config):
      ProditionsCache("StrawResponse",config.verbose()),
      _useDb(config.useDb()),_maker(config) {}

    void initialize() {
      _strawDrift_p = std::make_unique<ProditionsHandle<StrawDrift> >();
      _strawElectronics_p = std::make_unique<ProditionsHandle<StrawElectronics> >();
      _strawPhysics_p = std::make_unique<ProditionsHandle<StrawPhysics> >();
    }
    set_t makeSet(art::EventID const& eid) {
      auto sd = _strawDrift_p->get(eid);
      auto se = _strawElectronics_p->get(eid);
      auto sp = _strawPhysics_p->get(eid);
      auto ss = sd.getCids();
      ss.merge(set_t(se.getCids()));
      ss.merge(set_t(sp.getCids()));
      return ss;
    }
    DbIoV makeIov(art::EventID const& eid) {
      _strawDrift_p->get(eid);
      _strawElectronics_p->get(eid);
      _strawPhysics_p->get(eid);
      auto iov = _strawDrift_p->iov();
      iov.overlap(_strawElectronics_p->iov());
      iov.overlap(_strawPhysics_p->iov());
      return iov;
    }
    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      auto sd = _strawDrift_p->getPtr(eid);
      auto se = _strawElectronics_p->getPtr(eid);
      auto sp = _strawPhysics_p->getPtr(eid);
      return _maker.fromFcl(sd,se,sp);
    }

  private:
    bool _useDb;
    StrawResponseMaker _maker;

    // these handles are not default constructed
    // so to not create a dependency loop on construction
    std::unique_ptr<ProditionsHandle<StrawDrift> > _strawDrift_p;
    std::unique_ptr<ProditionsHandle<StrawElectronics> > _strawElectronics_p;
    std::unique_ptr<ProditionsHandle<StrawPhysics> > _strawPhysics_p;
  };
};

#endif
