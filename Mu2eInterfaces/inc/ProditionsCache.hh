#ifndef Mu2eInterfaces_ProditionsCache_hh
#define Mu2eInterfaces_ProditionsCache_hh
#include <memory>
#include <tuple>
#include <string>
#include <set>
#include <shared_mutex>
#include <mutex>
#include <chrono>
#include <iostream>

#include "canvas/Persistency/Provenance/EventID.h"
#include "Offline/DbTables/inc/DbIoV.hh"
#include "Offline/Mu2eInterfaces/inc/ProditionsEntity.hh"

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

    // what is actually held in the cache
    // new iovs are added if the data is the same
    struct cacheItem {
      ProditionsEntity::ptr _p;
      std::vector<DbIoV> _iovs;
    };

    ProditionsCache(std::string name, int verbose=0):
      _name(name),_verbose(verbose),_initialized(false) {}
    virtual ~ProditionsCache() {}

    // the following are provided by the
    // concrete class
    //virtual std::string const& name() const =0 ;
    std::string const& name() const { return _name;}
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
      bool found = false;
      ProditionsEntity::ptr p;
      set_t cids;
      DbIoV iov;
      { // start read lock scope
        std::shared_lock lock(_mutex);
        p = findByRun(eid,iov);  // if found, iov is valid
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
         p = findByRun(eid,iov);
         if(!p) {

           p = makeEntity(eid);
           cids = makeSet(eid);
           p->addCids(cids);
           iov = makeIov(eid);

           // at this point, we might have existing cache items with
           // the same cids, but not the relevant iov,
           // in this case just add the iov
           for(auto& ci : _cache) {
             if(ci._p->getCids()==cids) {
               ci._iovs.emplace_back(iov);
               found = true;
               break;
             }
           }
           if(!found) {
             cacheItem ci;
             ci._p = p;
             ci._iovs.emplace_back(iov);
             _cache.emplace_back(ci);
             made = true;
             if(_verbose>7) p->print(std::cout);
           }
         } // p not found

         auto etime = std::chrono::high_resolution_clock::now();
         dt = std::chrono::duration_cast<std::chrono::microseconds>
                                               ( etime - mtime );
        _lockTime += dt;  // time we spent write locked

      } // endif not in cache, write lock now destroyed

      if(_verbose>1) {
        if(made) {
          if(found) {
            std::cout<< "ProditionsCache::update made new iov for "
                     << name() << std::endl;
          } else {
            std::cout<< "ProditionsCache::update made new "
                     << name() << std::endl;
          }
        } else {
          std::cout<< "ProditionsCache::update return cached "<< name() << std::endl;
        }
        std::cout << "     iov " << iov.to_string(true);
        std::cout << "     cids ";
        for(auto cid : cids) std::cout << cid << " " ;
        std::cout << std::endl;
      }

      return std::make_tuple(p,iov);

    } // end update

    // is there a cache entry covering this run/subrun?
    // return good pointer or null, and fill iov
    ProditionsEntity::ptr  findByRun(art::EventID eid, DbIoV& iov) {
      uint32_t run = eid.run();
      uint32_t subrun = eid.subRun();
      for(auto const& ci : _cache) {
        for(auto const& ii : ci._iovs) {
          if(ii.inInterval(run,subrun)) {
            iov = ii;
            return ci._p;
          }
        }
      }
      return ProditionsEntity::ptr();
    }

  private:
    std::string _name;
    int _verbose;
    bool _initialized;
    std::vector<cacheItem> _cache;

  };

}

#endif
