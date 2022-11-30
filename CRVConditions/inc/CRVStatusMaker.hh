#ifndef CRVConditions_CRVStatusMaker_hh
#define CRVConditions_CRVStatusMaker_hh

//
// construct a CRVStatus proditions entity
// (flags for bad channels) from fcl or database
//

#include "Offline/CRVConditions/inc/CRVStatus.hh"
#include "Offline/CRVConfig/inc/CRVStatusConfig.hh"
#include "Offline/DbTables/inc/CRVBadChan.hh"

namespace mu2e {

class CRVStatusMaker {
 public:
  CRVStatusMaker(CRVStatusConfig const& config) : _config(config) {}

  CRVStatus::ptr_t fromFcl();
  CRVStatus::ptr_t fromDb(CRVBadChan::cptr_t cbc_p);

 private:
  // this object needs to be thread safe,
  // _config should only be initialized once
  const CRVStatusConfig _config;
};
}  // namespace mu2e

#endif
