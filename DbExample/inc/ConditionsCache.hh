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
    typedef ConditionsEntity2::set_t set_t;

    virtual ~ConditionsCache() {}
    virtual std::string const& name() const =0 ;
    virtual ret_t update(art::EventID const& eid) =0;

    void push(ConditionsEntity2::ptr const& p);
    ConditionsEntity2::ptr find(set_t const& s);

  private:
    std::vector<ConditionsEntity2::ptr> _cache;

  };

}

#endif
