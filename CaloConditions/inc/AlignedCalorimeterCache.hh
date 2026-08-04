#ifndef TrackerConditions_AlignedCalorimeterCache_hh
#define TrackerConditions_AlignedCalorimeterCache_hh
//
// hold a set of run-dependent conditions objects
// and update them when needed
//
#include "Offline/Mu2eInterfaces/inc/ProditionsCache.hh"
#include "Offline/DbService/inc/DbHandle.hh"
#include "Offline/DbTables/inc/CalAlignElement.hh"
#include "Offline/CaloConditions/inc/AlignedCalorimeterMaker.hh"


namespace mu2e {

  class AlignedCalorimeterCache : public ProditionsCache {
    public:
      AlignedCalorimeterCache(AlignedCalConfig const& config):
        ProditionsCache(Calorimeter::cxname,config.verbose()),
        _useDb(config.useDb()),_maker(config) {}

      void initialize() {
        if (_useDb) {
          _tadisk_p = std::make_unique<DbHandle<CalAlignDisk>>();
          _tacrys_p = std::make_unique<DbHandle<CalAlignCrystal>>();
        }
      }

      set_t makeSet(art::EventID const& eid) {
        // get the tables up to date
        ProditionsEntity::set_t cids;
        if (_useDb) {
          _tadisk_p->get(eid);
          _tacrys_p->get(eid);
          cids.insert(_tadisk_p->cid());
          cids.insert(_tacrys_p->cid());
        }
        return cids;
      }

      DbIoV makeIov(art::EventID const& eid) {
        DbIoV iov;
        iov.setMax(); // start with full IOV range
        if (_useDb) {
          iov.overlap(_tadisk_p->iov());
          iov.overlap(_tacrys_p->iov());
        }
        return iov;
      }

      ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
        if(_useDb) {
          return _maker.fromDb(_tadisk_p->getPtr(eid),_tacrys_p->getPtr(eid));
        } else {
          return _maker.fromFcl();
        }
      }


    private:
      bool _useDb;
      AlignedCalorimeterMaker _maker;

      std::unique_ptr<DbHandle<CalAlignDisk>>    _tadisk_p;
      std::unique_ptr<DbHandle<CalAlignCrystal>> _tacrys_p;
  };
}
#endif
