#ifndef AnalysisConditions_MVACatalogCache_hh
#define AnalysisConditions_MVACatalogCache_hh

//
// MVACatalogCache for ProditionsCache
// This is a templated class that takes the MVAStruct<DETAIL> class and the MVAToolDb class
// For an example use see AnalysisConditions/inc/TrkQualCatalogCache.hh
//

#include "Mu2eInterfaces/inc/ProditionsCache.hh"
#include "DbTables/inc/DbIoV.hh"
#include "DbService/inc/DbHandle.hh"
#include "DbTables/inc/MVAToolDb.hh"
#include "AnalysisConditions/inc/MVACatalogMaker.hh"

namespace mu2e {
  template <class T, class DB>
  class MVACatalogCache : public ProditionsCache {
  public: 
    MVACatalogCache(std::string name, MVACatalogConfig const& config):
      ProditionsCache(name,config.verbose()),
      _useDb(config.useDb()),_maker(config) {}

    void initialize() {
      if(_useDb) {
	_tqDb_p  = std::make_unique<DbHandle<DB> >();
      }
    }
    
    set_t makeSet(art::EventID const& eid) {
      ProditionsEntity::set_t cids;
      if(_useDb) { // use fcl config, overwrite part from DB
	// get the tables up to date
	_tqDb_p->get(eid);
	// save which data goes into this instance of the service
	cids.insert(_tqDb_p->cid());
      }
      return cids;
    }
    
    DbIoV makeIov(art::EventID const& eid) {
      DbIoV iov;
      iov.setMax(); // start with full IOV range
      if(_useDb) { // use fcl config, overwrite part from DB
	// get the tables up to date
	_tqDb_p->get(eid);
	// restrict the valid range ot the overlap
	iov.overlap(_tqDb_p->iov());
      }
      return iov;
    }
    
    ProditionsEntity::ptr makeEntity(art::EventID const& eid) {
      if(_useDb) {
	return _maker.fromDb( _tqDb_p->getPtr(eid) );
      } else {
	return _maker.fromFcl();
      }
    }
    
  private:
    bool _useDb;
    MVACatalogMaker<T, DB> _maker;

    // these handles are not default constructed
    // so the db can be completely turned off
    std::unique_ptr<DbHandle<DB> > _tqDb_p;
  };
};

#endif
