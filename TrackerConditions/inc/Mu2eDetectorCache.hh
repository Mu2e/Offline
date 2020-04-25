#ifndef TrackerConditions_Mu2eDetectorCache_hh
#define TrackerConditions_Mu2eDetectorCache_hh

//
// hold a set of run-dependent conditions objects
// and update them when needed
//

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "TrackerConditions/inc/Mu2eDetectorMaker.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"

namespace mu2e {
  class Mu2eDetectorCache : public ProditionsCache {
  public: 
    Mu2eDetectorCache(Mu2eDetectorConfig const& config):
      ProditionsCache("Mu2eDetector",config.verbose()),
      _useDb(config.useDb()),_maker(config) {}


    void initialize() {
      _mu2eMaterial_p = std::make_unique<ProditionsHandle<Mu2eMaterial> >();
      _alignedDetector_p = std::make_unique<ProditionsHandle<Tracker> >();
    }

    set_t makeSet(art::EventID const& eid) {
      auto mm = _mu2eMaterial_p->getPtr(eid);
      auto tr = _alignedDetector_p->getPtr(eid);
      // this is the set of DB cid's they depend on
      ProditionsEntity::set_t ss;
      ss.merge(set_t(mm->getCids()));
      ss.merge(set_t(tr->getCids()));
      return ss;
    }

    DbIoV makeIov(art::EventID const& eid) {
      _mu2eMaterial_p->get(eid); // check up to date
      _alignedDetector_p->get(eid);
      auto iov = _mu2eMaterial_p->iov();
      iov.overlap(_alignedDetector_p->iov());
      return iov;
    }

    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      auto mm = _mu2eMaterial_p->getPtr(eid);
      auto tr = _alignedDetector_p->getPtr(eid);
      return _maker.fromFcl(mm,tr);
    }

  private:
    bool _useDb;
    Mu2eDetectorMaker _maker;

    // this handle is not default constructed
    // so to not create a dependency loop on construction
    std::unique_ptr<ProditionsHandle<Mu2eMaterial> > _mu2eMaterial_p;
    std::unique_ptr<ProditionsHandle<Tracker> > _alignedDetector_p;

  };
};

#endif
