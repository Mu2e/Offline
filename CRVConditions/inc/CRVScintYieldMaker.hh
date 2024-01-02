#ifndef CRVConditions_CRVScintYieldMaker_hh
#define CRVConditions_CRVScintYieldMaker_hh

//
// construct a CRVScintYield proditions entity
// (scintillation yield spread of CRV counters)
// from fcl or database
//

#include "Offline/CRVConditions/inc/CRVScintYield.hh"
#include "Offline/CRVConfig/inc/CRVScintYieldConfig.hh"
#include "Offline/DbTables/inc/CRVScint.hh"

namespace mu2e {

class CRVScintYieldMaker {
 public:
  CRVScintYieldMaker(CRVScintYieldConfig const& config) : _config(config) {}

  CRVScintYield::ptr_t fromFcl();
  CRVScintYield::ptr_t fromDb(CRVScint::cptr_t sci_p);

 private:
  // this object needs to be thread safe,
  // _config should only be initialized once
  const CRVScintYieldConfig _config;
};
}  // namespace mu2e

#endif
