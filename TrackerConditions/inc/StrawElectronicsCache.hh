#ifndef TrackerConditions_StrawElectronicsCache_hh
#define TrackerConditions_StrawElectronicsCache_hh

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
#include "DbTables/inc/DbIoV.hh"
#include "DbService/inc/DbHandle.hh"
#include "DbTables/inc/TrkDelayPanel.hh"
#include "DbTables/inc/TrkPreampRStraw.hh"
#include "DbTables/inc/TrkPreampStraw.hh"
#include "DbTables/inc/TrkThresholdRStraw.hh"

#include "TrackerConditions/inc/StrawElectronicsMaker.hh"


namespace mu2e {
  class StrawElectronicsCache : public ProditionsCache {
  public: 
    StrawElectronicsCache(StrawElectronicsConfig const& config):
      _name("StrawElectronics"),_maker(config),
      _verbose(config.verbose()),_useDb(config.useDb()) {}

    std::string const& name() const { return _name; }

    ProditionsCache::ret_t update(art::EventID const& eid) {

      // lock access to the data, will release when this method returns
      LockGuard lock(*this);

      if(_verbose>1) std::cout << "in cache useDb =  "<< _useDb<< std::endl;
      // delayed construction here so that the handles
      // inside the service can reference the service 
      if(_useDb && !_tdp_p) {
	_tdp_p  = std::make_unique<DbHandle<TrkDelayPanel>>();
	_tprs_p = std::make_unique<DbHandle<TrkPreampRStraw>>();
	_tps_p  = std::make_unique<DbHandle<TrkPreampStraw>>();
	_ttrs_p = std::make_unique<DbHandle<TrkThresholdRStraw>>();
      }

      // figure out what versions of tables are needed
      ProditionsEntity::set_t cids;
      DbIoV iov;
      iov.setMax(); // start with full IOV range
      if(_useDb) { // use fcl config, overwrite part from DB
	// get the tables up to date
	_tdp_p->get(eid);
	_tprs_p->get(eid);
	_tps_p->get(eid);
	_ttrs_p->get(eid);
	// save which data goes into this instance of the service
	cids.insert(_tdp_p->cid());
	cids.insert(_tprs_p->cid());
	cids.insert(_tps_p->cid());
	cids.insert(_ttrs_p->cid());
	// restrict the valid range ot the overlap
	iov.overlap(_tdp_p->iov());
	iov.overlap(_tprs_p->iov());
	iov.overlap(_tps_p->iov());
	iov.overlap(_ttrs_p->iov());
      }
      // here we would add dependence and IOV restrictions from 
      // other services this depends on
      // cids.insert(c1.getCids().begin(),c1.getCids().end());

      // see if this combination of tables is in the cache
      auto p = find(cids);

      if(_verbose>1) {
	if(!p) {
	  std::cout<< "making new StrawElectronics from cids=";
	} else {
	  std::cout<< "found StrawElectronics in cache with cids=";
	}
	for(auto x : cids) std::cout << x << " ";
	std::cout << std::endl;
      }

      if(!p) {
	if(_useDb) {
	  p = _maker.fromDb( _tdp_p->getPtr(eid),
			     _tprs_p->getPtr(eid), 
			     _tps_p->getPtr(eid),
			     _ttrs_p->getPtr(eid) );
        } else {
	  p = _maker.fromFcl();
	}
        p->addCids(cids);
        push(p);
      }

      return std::make_tuple(p,iov);
}

  private:
    std::string _name;
    StrawElectronicsMaker _maker;
    int _verbose;
    bool _useDb;

    // these handles are not default constructed
    // so the db can be completely turned off
    std::unique_ptr<DbHandle<TrkDelayPanel>> _tdp_p;
    std::unique_ptr<DbHandle<TrkPreampRStraw>> _tprs_p;
    std::unique_ptr<DbHandle<TrkPreampStraw>> _tps_p;
    std::unique_ptr<DbHandle<TrkThresholdRStraw>> _ttrs_p;

  };
};

#endif
