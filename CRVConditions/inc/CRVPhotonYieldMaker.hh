#ifndef CRVConditions_CRVPhotonYieldMaker_hh
#define CRVConditions_CRVPhotonYieldMaker_hh

//
// construct a CRVPhotonYield proditions entity
// (photon yield spread of CRV channels)
// from fcl or database
//

#include "Offline/CRVConditions/inc/CRVPhotonYield.hh"
#include "Offline/CRVConfig/inc/CRVPhotonYieldConfig.hh"
#include "Offline/DbTables/inc/CRVPhoton.hh"

namespace mu2e {

class CRVPhotonYieldMaker {
 public:
  CRVPhotonYieldMaker(CRVPhotonYieldConfig const& config) : _config(config) {}

  CRVPhotonYield::ptr_t fromFcl();
  CRVPhotonYield::ptr_t fromDb(CRVPhoton::cptr_t sci_p);

 private:
  // this object needs to be thread safe,
  // _config should only be initialized once
  const CRVPhotonYieldConfig _config;
};
}  // namespace mu2e

#endif
