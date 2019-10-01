#ifndef Mu2eInterfaces_ProditionsCache_hh
#define Mu2eInterfaces_ProditionsCache_hh
#include <memory>
#include <tuple>
#include <string>
#include <set>
#include <shared_mutex>
#include <mutex>
#include <chrono>

#include "canvas/Persistency/Provenance/EventID.h"
#include "DbTables/inc/DbIoV.hh"
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"

namespace mu2e {
  class ProditionsCache {

  protected:

    // lock for threaded access
    std::shared_mutex _mutex;

    // count the time waiting and locked
    std::chrono::microseconds _lockWaitTime;
    std::chrono::microseconds _lockTime;


  public: 
    typedef std::shared_ptr<ProditionsCache> ptr;
    typedef std::tuple<ProditionsEntity::ptr,DbIoV> ret_t;
    typedef ProditionsEntity::set_t set_t;

    ProditionsCache(std::string name, int verbose=0):
      _name(name),_verbose(verbose),_initialized(false) {}
    virtual ~ProditionsCache() {}

    // the following are provided by the 
    // concrete class
    //virtual std::string const& name() const =0 ;
    std::string const& name() { return _name;}
    // create pointers to services on demand
    // this allows lazy intialization and prevents
    // service dependence loops
    virtual void initialize() =0;
    // create the list of cid's of database cache objects we need
    virtual set_t makeSet(art::EventID const& eid) =0;
    // create the interval of validity
    virtual DbIoV makeIov(art::EventID const& eid) =0;
    // make a new entity, the data object itself
    virtual ProditionsEntity::ptr makeEntity(art::EventID const& eid) =0;

    // this is the main call to the cache asking for an existing
    // entity, creating and cacheing a new entity as needed
    ret_t update(art::EventID const& eid) {
      // do lazy initialization, don't bother with 
      // read lock since a bool can't be partially constructed
      if(!_initialized) {
	//gain write lock
	auto stime = std::chrono::high_resolution_clock::now();
	std::unique_lock lock(_mutex); // write lock
	auto mtime = std::chrono::high_resolution_clock::now();
	auto dt = std::chrono::duration_cast<std::chrono::microseconds>
                                               ( mtime - stime );
	 _lockWaitTime += dt;
	// check if another thread initialized while we were
	// waiting for write lock
	if(!_initialized) {
	  // derived class creates database and service dependencies
	  initialize(); 
	  _initialized = true;
	}
	auto etime = std::chrono::high_resolution_clock::now();
	dt = std::chrono::duration_cast<std::chrono::microseconds>
                                               ( etime - mtime );
	_lockTime += dt;  // time we spent write locked
      } // end initialize, write lock out of scope, released

      // gain shared read lock to find what 
      // set of tables are needed
      bool made = false;
      set_t cids;
      ProditionsEntity::ptr p;
      DbIoV iov;
      { // start read lock scope
	std::shared_lock lock(_mutex);
	// get the set of nubers that identifies the data
	cids = makeSet(eid);
	// look for it in the cache
	p = find(cids);
	if(p) { // also grab the iov while under this read lock
	  iov = makeIov(eid);
	}
      } // end read lock lifetime
      
      // if it was not found in cache, make it
      if(!p) {
	//gain write lock
	auto stime = std::chrono::high_resolution_clock::now();
	std::unique_lock lock(_mutex); // write lock
	auto mtime = std::chrono::high_resolution_clock::now();
	auto dt = std::chrono::duration_cast<std::chrono::microseconds>
                                               ( mtime - stime );
	 _lockWaitTime += dt;
	 // need to check again in case another thread made it
	 // between read lock and write lock
	 p = find(cids);
	 if(!p) {
	   p = makeEntity(eid); // make the data entity
	   p->addCids(cids); // label it
	   push(p); // put in the cache
	   made = true;
	   if(_verbose>2) p->print(std::cout);
	 }
	 iov = makeIov(eid); // new or old, iov is now valid
	 auto etime = std::chrono::high_resolution_clock::now();
	 dt = std::chrono::duration_cast<std::chrono::microseconds>
                                               ( etime - mtime );
	_lockTime += dt;  // time we spent write locked

      } // endif not in cache, write lock now destroyed
      
      if(_verbose>1) {
	if(made) {
	  std::cout<< "ProditionsCache::update made new "<< name() << std::endl;
	} else {
	  std::cout<< "ProditionsCache::update return cached "<< name() << std::endl;
	}
      }

      return std::make_tuple(p,iov);

    } // end update

    // put this object, with dependent set of CID's, in the cache
    void push(ProditionsEntity::ptr const& p) {
      _cache.emplace_back(p);
    }

    // is the object, with this set of CID's, 
    // which uniquely identifies it, in the cache?
    ProditionsEntity::ptr  find(set_t const& s) {
      for(auto const& ii : _cache) {
	if(ii->getCids()==s) return ii;
      }
      return ProditionsEntity::ptr();
    }
    
  private:
    std::string _name;
    int _verbose;
    bool _initialized;
    std::vector<ProditionsEntity::ptr> _cache;

  };

}

#endif
