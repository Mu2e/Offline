#ifndef ProditionService_ProditionsCache_hh
#define ProditionService_ProditionsCache_hh
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
      LockGuard(ProditionsCache& cache);
      ~LockGuard();
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

    void push(ProditionsEntity::ptr const& p);
    ProditionsEntity::ptr find(set_t const& s);

  private:
    std::vector<ProditionsEntity::ptr> _cache;

  };

}

#endif
