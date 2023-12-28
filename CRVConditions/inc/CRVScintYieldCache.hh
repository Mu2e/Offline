#ifndef CRVConditions_CRVScintYieldCache_hh
#define CRVConditions_CRVScintYieldCache_hh

#include "Offline/CRVConditions/inc/CRVScintYieldMaker.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsCache.hh"

namespace mu2e {

class CRVScintYieldCache : public ProditionsCache {
 public:
  CRVScintYieldCache(CRVScintYieldConfig const& config) :
      ProditionsCache(CRVScintYield::cxname, config.verbose()),
      _useDb(config.useDb()), _maker(config) {}

  void initialize() {}

  set_t makeSet(art::EventID const& eid) {
    ProditionsEntity::set_t cids;
    return cids;
  }

  DbIoV makeIov(art::EventID const& eid) {
    DbIoV iov;
    iov.setMax();
    return iov;
  }

  ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
    if (_useDb) {
      return _maker.fromDb();
    } else {
      return _maker.fromFcl();
    }
  }

 private:
  bool _useDb;
  CRVScintYieldMaker _maker;
};

}  // namespace mu2e

#endif
