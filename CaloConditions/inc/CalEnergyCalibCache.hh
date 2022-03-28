#ifndef TrackerConditions_CalEnergyCalibCache_hh
#define TrackerConditions_CalEnergyCalibCache_hh

#include "Offline/Mu2eInterfaces/inc/ProditionsCache.hh"
#include "Offline/CaloConditions/inc/CalEnergyCalibMaker.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"


namespace mu2e {
  class CalEnergyCalibCache : public ProditionsCache {
  public: 
    CalEnergyCalibCache(CalEnergyCalibConfig const& config):
      ProditionsCache(CalEnergyCalib::cxname,config.verbose()),
      _useDb(config.useDb()),_maker(config) {}

    void initialize() {
      _calenergycalib_p = std::make_unique<ProditionsHandle<CalEnergyCalib> >();
    }
    
    set_t makeSet(art::EventID const& eid) {
      ProditionsEntity::set_t cids;
      auto cal = _calenergycalib_p->get(eid);
      cids.insert(_calenergycalib_p->cid());
      return cids;
    }
    
    DbIoV makeIov(art::EventID const& eid) {
      _calenergycalib_p->get(eid);
      iov.overlap(_calenergycalib_p->iov());
      return iov;
    }
    
    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      auto cal = _calenergycalib_p->getPtr(eid);
      return _maker.fromFcl(cal);
    }

  private:
    bool _useDb;
    CalEnergyCalibMaker _maker;
    std::unique_ptr<ProditionsHandle<CalEnergyCalib> > _calenergycalib_p;
  };
};

#endif
