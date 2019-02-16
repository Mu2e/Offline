#ifndef Mu2eInterfaces_ProditionsCache_hh
#define Mu2eInterfaces_ProditionsCache_hh
#include <memory>
#include <tuple>
#include <string>
#include <set>
#include <mutex>
#include <chrono>

#include "canvas/Persistency/Provenance/EventID.h"
#include "DbTables/inc/DbIoV.hh"
#include "Mu2eInterfaces/inc/ProditionsEntity.hh"

namespace mu2e {
  class ProditionsCache {

  protected:
    // for subclass to make sure that when the update 
    // method exits, the lock is released
    class LockGuard {
      ProditionsCache& _cache;
    public:
      LockGuard(ProditionsCache& cache):_cache(cache) {
	_cache._lock.lock();
	_cache._lockTime = std::chrono::high_resolution_clock::now();
      }

      ~LockGuard() {
	auto end_time = std::chrono::high_resolution_clock::now();
	auto delta = std::chrono::duration_cast<std::chrono::microseconds>
	  (end_time - _cache._lockTime);
	_cache._lockTotalTime += delta;
	_cache._lock.unlock();
      }
    };

    // lock for threaded access
    std::mutex _lock;

    // count the time locked
    std::chrono::high_resolution_clock::time_point _lockTime;
    std::chrono::microseconds _lockTotalTime;


  public: 
    typedef std::shared_ptr<ProditionsCache> ptr;
    typedef std::tuple<ProditionsEntity::ptr,DbIoV> ret_t;
    typedef ProditionsEntity::set_t set_t;

    virtual ~ProditionsCache() {}
    virtual std::string const& name() const =0 ;
    virtual ret_t update(art::EventID const& eid) =0;

    // put this object, with dependent set of CID's, in the cache
    void push(ProditionsEntity::ptr const& p) {
      _cache.emplace_back(p);
    }

    // is the object, with this set of CID's, in the cache?
    ProditionsEntity::ptr  find(set_t const& s) {
      for(auto const& ii : _cache) {
	if(ii->getCids()==s) return ii;
      }
      return ProditionsEntity::ptr();
    }
    
  private:
    std::vector<ProditionsEntity::ptr> _cache;

  };

}

#endif
