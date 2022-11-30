#ifndef CRVConditions_CRVCalibMaker_hh
#define CRVConditions_CRVCalibMaker_hh

//
// construct a CRVCalib proditions entity (CRV SiPM pedestal, energy, time)
// from fcl or database
//

#include "Offline/CRVConditions/inc/CRVCalib.hh"
#include "Offline/CRVConfig/inc/CRVCalibConfig.hh"
#include "Offline/DbTables/inc/CRVSiPM.hh"
#include "Offline/DbTables/inc/CRVTime.hh"

namespace mu2e {

class CRVCalibMaker {
 public:
  CRVCalibMaker(CRVCalibConfig const& config) : _config(config) {}

  CRVCalib::ptr_t fromFcl();
  CRVCalib::ptr_t fromDb(CRVSiPM::cptr_t sip_p, CRVTime::cptr_t tim_p);

 private:
  // this object needs to be thread safe,
  // _config should only be initialized once
  const CRVCalibConfig _config;
};
}  // namespace mu2e

#endif
