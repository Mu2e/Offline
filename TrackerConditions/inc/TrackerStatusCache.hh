#ifndef TrackerConditions_TrackerStatusCache_hh
#define TrackerConditions_TrackerStatusCache_hh

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
#include "DbService/inc/DbHandle.hh"
#include "DbTables/inc/TrkElementStatus.hh"
#include "TrackerConditions/inc/TrackerStatusMaker.hh"


namespace mu2e {
  class TrackerStatusCache : public ProditionsCache {
  public: 
    TrackerStatusCache(TrackerStatusConfig const& config):
      ProditionsCache(TrackerStatus::cxname,config.settings().verbose()),
      _useDb(config.settings().useDb()),_maker(config) {}


    void initialize() {
      if(_useDb) {
	_tpls_p = std::make_unique<DbHandle<TrkPlaneStatus>>();
	_tpas_p = std::make_unique<DbHandle<TrkPanelStatus>>();
	_tssl_p = std::make_unique<DbHandle<TrkStrawStatusLong>>(); 
	_tsss_p = std::make_unique<DbHandle<TrkStrawStatusShort>>();
      }
    }

    set_t makeSet(art::EventID const& eid) {
      // get the tables up to date
      ProditionsEntity::set_t cids;
      if(_useDb) {
	_tpls_p->get(eid);
	_tpas_p->get(eid);
	_tssl_p->get(eid);
	_tsss_p->get(eid);
	cids.insert(_tpls_p->cid());
	cids.insert(_tpas_p->cid());
	cids.insert(_tssl_p->cid());
	cids.insert(_tsss_p->cid());
      }
      return cids;
    }

    DbIoV makeIov(art::EventID const& eid) {
      DbIoV iov;
      iov.setMax(); // start with full IOV range
      if(_useDb) {
	iov.overlap(_tpls_p->iov());
	iov.overlap(_tpas_p->iov());
	iov.overlap(_tssl_p->iov());
	iov.overlap(_tsss_p->iov());
      }
      return iov;
    }

    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      if(_useDb) {
	return _maker.fromDb( _tpls_p->getPtr(eid),
	    _tpas_p->getPtr(eid), 
	    _tssl_p->getPtr(eid),
	    _tsss_p->getPtr(eid) );
      } else {
	return _maker.fromFcl();
      }
    }


  private:
    bool _useDb;
    TrackerStatusMaker _maker;

    std::unique_ptr<DbHandle<TrkPlaneStatus>>       _tpls_p;
    std::unique_ptr<DbHandle<TrkPanelStatus>>       _tpas_p;
    std::unique_ptr<DbHandle<TrkStrawStatusLong>>   _tssl_p;
    std::unique_ptr<DbHandle<TrkStrawStatusShort>>  _tsss_p;
  };
};

#endif
