#ifndef DbExample_ConditionsCache_hh
#define DbExample_ConditionsCache_hh
#include <memory>
#include <tuple>
#include <string>
#include <set>
#include <mutex>
#include <chrono>

#include "canvas/Persistency/Provenance/EventID.h"
#include "DbTables/inc/DbIoV.hh"
#include "DbExample/inc/ConditionsEntity2.hh"

namespace mu2e {
  class ConditionsCache {

  protected:
    // for subclass to make sure that when the update 
    // method exits, the lock is released
    class LockGuard {
      ConditionsCache& _cache;
    public:
      LockGuard(ConditionsCache& cache);
      ~LockGuard();
    };
    // lock for threaded access
    std::mutex _lock;
    // count the time locked
    std::chrono::high_resolution_clock::time_point _lockTime;
    std::chrono::microseconds _lockTotalTime;


  public: 
    typedef std::shared_ptr<ConditionsCache> ptr;
    typedef std::tuple<ConditionsEntity2::ptr,DbIoV> ret_t;
    typedef std::set<int> set_t;

    class item {
    public:
      item(ConditionsEntity2::ptr p, set_t s):_p(p),_s(s) {}
      ConditionsEntity2::ptr const& ptr() const { return _p; }
      set_t const& set() const { return _s; }
    private:
      ConditionsEntity2::ptr _p;
      set_t _s;
    };

    virtual ~ConditionsCache() {}
    virtual std::string const& name() const =0 ;
    virtual ret_t update(art::EventID const& eid) =0;

    void push(ConditionsEntity2::ptr const& p, set_t const& s);
    ConditionsEntity2::ptr find(set_t const& s);

    std::vector<item> _cache;

  };

}

#endif
