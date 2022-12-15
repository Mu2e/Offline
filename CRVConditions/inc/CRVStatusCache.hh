#ifndef CRVConditions_CRVStatusCache_hh
#define CRVConditions_CRVStatusCache_hh

#include "Offline/CRVConditions/inc/CRVStatusMaker.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsCache.hh"

namespace mu2e {

class CRVStatusCache : public ProditionsCache {
 public:
  CRVStatusCache(CRVStatusConfig const& config) :
      ProditionsCache(CRVStatus::cxname, config.verbose()),
      _useDb(config.useDb()), _maker(config) {}

  void initialize() {
    if (_useDb) {
      _cbc_p = std::make_unique<DbHandle<CRVBadChan>>();
    }
  }

  set_t makeSet(art::EventID const& eid) {
    ProditionsEntity::set_t cids;
    if (_useDb) {
      _cbc_p->get(eid);
      cids.insert(_cbc_p->cid());
    }
    return cids;
  }

  DbIoV makeIov(art::EventID const& eid) {
    DbIoV iov;
    iov.setMax();
    if (_useDb) {
      _cbc_p->get(eid);
      iov.overlap(_cbc_p->iov());
    }
    return iov;
  }

  ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
    if (_useDb) {
      return _maker.fromDb(_cbc_p->getPtr(eid));
    } else {
      return _maker.fromFcl();
    }
  }

 private:
  bool _useDb;
  CRVStatusMaker _maker;

  // these handles are not default constructed
  // so the db can be completely turned off
  std::unique_ptr<DbHandle<CRVBadChan>> _cbc_p;
};

}  // namespace mu2e

#endif
