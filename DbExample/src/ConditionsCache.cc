#include "DbExample/inc/ConditionsCache.hh"

namespace mu2e {

  void ConditionsCache::push(ConditionsEntity2::ptr const& p, set_t const& s) {
    _cache.emplace_back(p,s);
  }
  ConditionsEntity2::ptr ConditionsCache::find(set_t const& s) {
    for(auto const& ii : _cache) {
      if(ii.set()==s) {
	return ii.ptr();
      }
    }
    return ConditionsEntity2::ptr();
  }

}
