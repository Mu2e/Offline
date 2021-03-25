#ifndef CaloConditions_CaloDAQConditionsCache_hh
#define CaloConditions_CaloDAQConditionsCache_hh

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
//#include "DbTables/inc/DbIoV.hh"
#include "DbService/inc/DbHandle.hh"
#include "DbTables/inc/DIRACtoCalo.hh"
#include "DbTables/inc/CalotoDIRAC.hh"

#include "CaloConditions/inc/CaloDAQConditionsMaker.hh"


namespace mu2e {
  class CaloDAQConditionsCache : public ProditionsCache {
  public: 
    CaloDAQConditionsCache(CaloDAQConditionsConfig const& config):
    ProditionsCache(CaloDAQConditions::cxname,config.verbose()),
      _useDb(config.useDb()),_maker(config) {}

    void initialize() {
      if(_useDb) {
	_tdtc_p = std::make_unique<DbHandle<DIRACtoCalo>>();
	_tctd_p = std::make_unique<DbHandle<CalotoDIRAC>>();
      }
    }
    
    set_t makeSet(art::EventID const& eid) {
      ProditionsEntity::set_t cids;
      if(_useDb) { // use fcl config, overwrite part from DB
	// get the tables up to date
	_tdtc_p->get(eid);
	_tctd_p->get(eid);
	// save which data goes into this instance of the service
	cids.insert(_tdtc_p->cid());
	cids.insert(_tctd_p->cid());
      }
      return cids;
    }
    
    DbIoV makeIov(art::EventID const& eid) {
      DbIoV iov;
      iov.setMax(); // start with full IOV range
      if(_useDb) { // use fcl config, overwrite part from DB
	// get the tables up to date
	_tdtc_p->get(eid);
	_tctd_p->get(eid);
	// restrict the valid range ot the overlap
	iov.overlap(_tdtc_p->iov());
	iov.overlap(_tctd_p->iov());
      }
      return iov;
    }
    
    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      if(_useDb) {
	return _maker.fromDb( _tdtc_p->getPtr(eid),
			      _tctd_p->getPtr(eid) ); 
      } else {
	return _maker.fromFcl();
      }
    }
    
  private:
    bool _useDb;
    CaloDAQConditionsMaker _maker;

    // these handles are not default constructed
    // so the db can be completely turned off
    std::unique_ptr<DbHandle<DIRACtoCalo>> _tdtc_p;
    std::unique_ptr<DbHandle<CalotoDIRAC>> _tctd_p;

  };
};

#endif
