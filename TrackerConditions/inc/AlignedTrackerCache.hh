#ifndef TrackerConditions_AlignedTrackerCache_hh
#define TrackerConditions_AlignedTrackerCache_hh

//
// hold a set of run-dependent conditions objects
// and update them when needed
//

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
#include "DbService/inc/DbHandle.hh"
#include "DbTables/inc/TrkAlignTracker.hh"
#include "DbTables/inc/TrkAlignPlane.hh"
#include "DbTables/inc/TrkAlignPanel.hh"
#include "TrackerConditions/inc/AlignedTrackerMaker.hh"


namespace mu2e {
  class AlignedTrackerCache : public ProditionsCache {
  public: 
    AlignedTrackerCache(AlignedTrackerConfig const& config):
      _name("AlignedTracker"),_maker(config),
      _verbose(config.verbose()),_useDb(config.useDb()) {}

    std::string const& name() const { return _name; }

    ProditionsCache::ret_t update(art::EventID const& eid) {

      // lock access to the data, will release when this method returns
      LockGuard lock(*this);

      if(_verbose>1) std::cout << "AlignedTracker useDb =  "
			       << _useDb<< std::endl;
      // delayed construction here so that the handles
      // inside the service can reference the service 
      if(_useDb && !_tatr_p) {
	_tatr_p = std::make_unique<DbHandle<TrkAlignTracker>>();
	_tapl_p = std::make_unique<DbHandle<TrkAlignPlane>>();
	_tapa_p = std::make_unique<DbHandle<TrkAlignPanel>>();
      }

      // figure out what versions of tables are needed
      ProditionsEntity::set_t cids;
      DbIoV iov;
      iov.setMax(); // start with full IOV range
      if(_useDb) { // use fcl config, overwrite part from DB
	// get the tables up to date
	_tatr_p->get(eid);
	_tapl_p->get(eid);
	_tapa_p->get(eid);
	// save which data goes into this instance of the service
	cids.insert(_tatr_p->cid());
	cids.insert(_tapl_p->cid());
	cids.insert(_tapa_p->cid());
	// restrict the valid range to the overlap
	iov.overlap(_tatr_p->iov());
	iov.overlap(_tapl_p->iov());
	iov.overlap(_tapa_p->iov());
      }
      // here we would add dependence and IOV restrictions from 
      // other services this depends on
      // cids.insert(c1.getCids().begin(),c1.getCids().end());

      // see if this combination of tables is in the cache
      auto p = find(cids);

      if(_verbose>1) {
	if(!p) {
	  std::cout<< "making new AlignedTracker from cids=";
	} else {
	  std::cout<< "found AlignedTracker in cache with cids=";
	}
	for(auto x : cids) std::cout << x << " ";
	std::cout << std::endl;
      }

      if(!p) {
	if(_useDb) {
	  p = _maker.fromDb( _tatr_p->getPtr(eid),
			     _tapl_p->getPtr(eid), 
			     _tapa_p->getPtr(eid) );
        } else {
	  p = _maker.fromFcl();
	}
        p->addCids(cids);
        push(p);

	if(_verbose>2) p->print(std::cout);

      }

      return std::make_tuple(p,iov);
    }

  private:
    std::string _name;
    AlignedTrackerMaker _maker;
    int _verbose;
    bool _useDb;

    std::unique_ptr<DbHandle<TrkAlignTracker>> _tatr_p;
    std::unique_ptr<DbHandle<TrkAlignPlane>>   _tapl_p;
    std::unique_ptr<DbHandle<TrkAlignPanel>>   _tapa_p;

  };
};

#endif
