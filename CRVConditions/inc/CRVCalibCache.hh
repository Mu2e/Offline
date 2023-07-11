#ifndef CRVConditions_CRVCalibCache_hh
#define CRVConditions_CRVCalibCache_hh

#include "Offline/CRVConditions/inc/CRVCalibMaker.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsCache.hh"

namespace mu2e {

class CRVCalibCache : public ProditionsCache {
 public:
  CRVCalibCache(CRVCalibConfig const& config) :
      ProditionsCache(CRVCalib::cxname, config.verbose()),
      _useDb(config.useDb()), _maker(config) {}

  void initialize() {
    if (_useDb) {
      _sip_p = std::make_unique<DbHandle<CRVSiPM>>();
      _tim_p = std::make_unique<DbHandle<CRVTime>>();
    }
  }

  set_t makeSet(art::EventID const& eid) {
    ProditionsEntity::set_t cids;
    if (_useDb) {
      _sip_p->get(eid);
      _tim_p->get(eid);
      cids.insert(_sip_p->cid());
      cids.insert(_tim_p->cid());
    }
    return cids;
  }

  DbIoV makeIov(art::EventID const& eid) {
    DbIoV iov;
    iov.setMax();
    if (_useDb) {
      _sip_p->get(eid);
      _tim_p->get(eid);
      iov.overlap(_sip_p->iov());
      iov.overlap(_tim_p->iov());
    }
    return iov;
  }

  ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
    if (_useDb) {
      return _maker.fromDb(_sip_p->getPtr(eid), _tim_p->getPtr(eid));
    } else {
      return _maker.fromFcl();
    }
  }

 private:
  bool _useDb;
  CRVCalibMaker _maker;

  // these handles are not default constructed
  // so the db can be completely turned off
  std::unique_ptr<DbHandle<CRVSiPM>> _sip_p;
  std::unique_ptr<DbHandle<CRVTime>> _tim_p;
};

}  // namespace mu2e

#endif
