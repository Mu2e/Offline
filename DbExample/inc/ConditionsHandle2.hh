#ifndef DbExample_ConditionsHandle2_hh
#define DbExample_ConditionsHandle2_hh

//
// A safe pointer to a ConditionsEntity2
//

#include <string>
#include "canvas/Persistency/Provenance/EventID.h"
#include "DbExample/inc/ConditionsService2.hh"
#include "DbTables/inc/DbIoV.hh"

namespace mu2e {
  template <typename ENTITY>
  class ConditionsHandle2
  {
  public:
    ConditionsHandle2() {
      art::ServiceHandle<ConditionsService2> sg;
      ENTITY e;
      _cptr = sg->getCache(e.name());
    }
    ~ConditionsHandle2() { }

    ENTITY const& get(art::EventID const& eid) { 
      uint32_t r = eid.run();
      uint32_t s = eid.subRun();
      if(!iov.inInterval(r,s)) {
	ConditionsEntity2::ptr bptr;
	std::tie(bptr,iov) = _cptr->update(eid);
	ptr = std::dynamic_pointer_cast
	  <const ENTITY,const ConditionsEntity2>(bptr);
      }
      return *ptr; 

    }

  private:
    ConditionsCache::ptr _cptr;
    std::shared_ptr<const ENTITY> ptr;
    DbIoV iov;
  };
}

#endif
