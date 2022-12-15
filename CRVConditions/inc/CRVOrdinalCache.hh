#ifndef CRVConditions_CRVOrdinalCache_hh
#define CRVConditions_CRVOrdinalCache_hh

#include "Offline/CRVConditions/inc/CRVOrdinalMaker.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsCache.hh"

namespace mu2e {

class CRVOrdinalCache : public ProditionsCache {
 public:
  CRVOrdinalCache(CRVOrdinalConfig const& config) :
      ProditionsCache(CRVOrdinal::cxname, config.verbose()),
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
  CRVOrdinalMaker _maker;
};

}  // namespace mu2e

#endif
