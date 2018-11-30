#include "DbExample/inc/ConditionsCache.hh"

namespace mu2e {

  void ConditionsCache::push(ConditionsEntity2::ptr const& p) {
    _cache.emplace_back(p);
  }
  ConditionsEntity2::ptr 
  ConditionsCache::find(set_t const& s) {
    for(auto const& ii : _cache) {
      if(ii->getCids()==s) return ii;
    }
    return ConditionsEntity2::ptr();
  }


  ConditionsCache::LockGuard::LockGuard(ConditionsCache& cache):_cache(cache) {
    _cache._lock.lock();
    _cache._lockTime = std::chrono::high_resolution_clock::now();
  }

  ConditionsCache::LockGuard::~LockGuard() {
    auto end_time = std::chrono::high_resolution_clock::now();
    auto delta = std::chrono::duration_cast<std::chrono::microseconds>
      (end_time - _cache._lockTime);
    _cache._lockTotalTime += delta;
    _cache._lock.unlock();
  }


}
