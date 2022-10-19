#ifndef STMConditions_STMEnergyCalibCache_hh
#define STMConditions_STMEnergyCalibCache_hh

#include "Offline/Mu2eInterfaces/inc/ProditionsCache.hh"
#include "Offline/STMConditions/inc/STMEnergyCalibMaker.hh"

namespace mu2e {

class STMEnergyCalibCache : public ProditionsCache {
 public:
  STMEnergyCalibCache(STMEnergyCalibConfig const& config) :
      ProditionsCache(STMEnergyCalib::cxname, config.verbose()),
      _useDb(config.useDb()), _maker(config) {}

  void initialize() {
    if (_useDb) {
      _sep_p = std::make_unique<DbHandle<STMEnergyPar>>();
      _ped_p = std::make_unique<DbHandle<STMPedestals>>();
    }
  }

  set_t makeSet(art::EventID const& eid) {
    ProditionsEntity::set_t cids;
    if (_useDb) {
      _sep_p->get(eid);
      _ped_p->get(eid);
      cids.insert(_sep_p->cid());
      cids.insert(_ped_p->cid());
    }
    return cids;
  }

  DbIoV makeIov(art::EventID const& eid) {
    DbIoV iov;
    iov.setMax();
    if (_useDb) {
      _sep_p->get(eid);
      _ped_p->get(eid);
      iov.overlap(_sep_p->iov());
      iov.overlap(_ped_p->iov());
    }
    return iov;
  }

  ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
    if (_useDb) {
      return _maker.fromDb(_sep_p->getPtr(eid), _ped_p->getPtr(eid));
    } else {
      return _maker.fromFcl();
    }
  }

 private:
  bool _useDb;
  STMEnergyCalibMaker _maker;

  // these handles are not default constructed
  // so the db can be completely turned off
  std::unique_ptr<DbHandle<STMEnergyPar>> _sep_p;
  std::unique_ptr<DbHandle<STMPedestals>> _ped_p;
};

}  // namespace mu2e

#endif
