#ifndef CaloConditions_CalCalibCache_hh
#define CaloConditions_CalCalibCache_hh


// This Proditions entitiy cache is for the combined calorimeter calibration output
// author: S. Middleton 2022

#include "Offline/Mu2eInterfaces/inc/ProditionsCache.hh"
#include "Offline/CaloConditions/inc/CalCalibMaker.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"


namespace mu2e {
  class CalCalibCache : public ProditionsCache {
  public:
    CalCalibCache(CalCalibConfig const& config):
      ProditionsCache(CalCalib::cxname,config.verbose()),
      _useDb(config.useDb()),_maker(config) {}

    void initialize() {
     if(_useDb) {
        _calenergycalib_p = std::make_unique<DbHandle<CalEnergyCalib>>();
        _caltimecalib_p = std::make_unique<DbHandle<CalTimeCalib>>();
      }
    }

    set_t makeSet(art::EventID const& eid) {
      ProditionsEntity::set_t cids;
      if(_useDb) {
        _calenergycalib_p->get(eid);
        _caltimecalib_p->get(eid);
        cids.insert(_calenergycalib_p->cid());
        cids.insert(_caltimecalib_p->cid());
      }
      return cids;
    }

    DbIoV makeIov(art::EventID const& eid) {
      DbIoV iov;
      iov.setMax(); // start with full IOV range

      if(_useDb) {
        _calenergycalib_p->get(eid);
        _caltimecalib_p->get(eid);

         iov.overlap(_calenergycalib_p->iov());
        iov.overlap(_caltimecalib_p->iov());
      }
      return iov;
    }

    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      if(_useDb) {
        return _maker.fromDb( _calenergycalib_p->get(eid),
                              _caltimecalib_p->get(eid));
      } else {
        return _maker.fromFcl();
      }
    }

  private:
    bool _useDb;
    CalCalibMaker _maker;
    std::unique_ptr<DbHandle<CalEnergyCalib>> _calenergycalib_p;
    std::unique_ptr<DbHandle<CalTimeCalib>> _caltimecalib_p;
  };
}

#endif
