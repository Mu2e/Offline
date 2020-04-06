#ifndef TrackerConditions_TrackerDAQConditionsCache_hh
#define TrackerConditions_TrackerDAQConditionsCache_hh

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
#include "DbTables/inc/DbIoV.hh"
#include "DbService/inc/DbHandle.hh"
#include "DbTables/inc/TrkDelayPanel.hh"
#include "DbTables/inc/TrkPreampRStraw.hh"
#include "DbTables/inc/TrkPreampStraw.hh"
#include "DbTables/inc/TrkThresholdRStraw.hh"

#include "TrackerConditions/inc/TrackerDAQConditionsMaker.hh"


namespace mu2e {
  class TrackerDAQConditionsCache : public ProditionsCache {
  public: 
    TrackerDAQConditionsCache(TrackerDAQConditionsConfig const& config):
      ProditionsCache("TrackerDAQConditions",config.verbose()),
      _useDb(config.useDb()),_maker(config) {}

    void initialize() {
      if(_useDb) {
	_tdts_p = std::make_unique<DbHandle<TrkDRACtoStraw>>();
	_trtp_p = std::make_unique<DbHandle<TrkROCtoPanel>>();
      }
    }
    
    set_t makeSet(art::EventID const& eid) {
      ProditionsEntity::set_t cids;
      if(_useDb) { // use fcl config, overwrite part from DB
	// get the tables up to date
	_tdts_p->get(eid);
	_trtp_p->get(eid);
	// save which data goes into this instance of the service
	cids.insert(_tdts_p->cid());
	cids.insert(_trtp_p->cid());
      }
      return cids;
    }
    
    DbIoV makeIov(art::EventID const& eid) {
      DbIoV iov;
      iov.setMax(); // start with full IOV range
      if(_useDb) { // use fcl config, overwrite part from DB
	// get the tables up to date
	_tdts_p->get(eid);
	_trtp_p->get(eid);
	// restrict the valid range ot the overlap
	iov.overlap(_tdts_p->iov());
	iov.overlap(_trtp_p->iov());
      }
      return iov;
    }
    
    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      if(_useDb) {
	return _maker.fromDb( _tdts_p->getPtr(eid),
			      _trtp_p->getPtr(eid) ); 
      } else {
	return _maker.fromFcl();
      }
    }
    
  private:
    bool _useDb;
    TrackerDAQConditionsMaker _maker;

    // these handles are not default constructed
    // so the db can be completely turned off
    std::unique_ptr<DbHandle<TrkDRACtoStraw>> _tdts_p;
    std::unique_ptr<DbHandle<TrkROCtoPanel>> _trtp_p;

  };
};

#endif
