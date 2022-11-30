#ifndef CRVConditions_CRVOrdinalMaker_hh
#define CRVConditions_CRVOrdinalMaker_hh

//
// construct a CRVOrdinal proditions entity
// (convert between CRV online and offline numbering)
// from fcl or database
//

#include "Offline/CRVConditions/inc/CRVOrdinal.hh"
#include "Offline/CRVConfig/inc/CRVOrdinalConfig.hh"

namespace mu2e {

class CRVOrdinalMaker {
 public:
  CRVOrdinalMaker(CRVOrdinalConfig const& config) : _config(config) {}

  CRVOrdinal::ptr_t fromFcl();
  CRVOrdinal::ptr_t fromDb();

 private:
  // this object needs to be thread safe,
  // _config should only be initialized once
  const CRVOrdinalConfig _config;
};
}  // namespace mu2e

#endif
