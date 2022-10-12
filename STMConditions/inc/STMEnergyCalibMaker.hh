#ifndef STMConditions_STMEnergyCalibMaker_hh
#define STMConditions_STMEnergyCalibMaker_hh

//
// construct a STMEnergyCalib proditions entity
// from fcl or database
//

#include "Offline/DbTables/inc/STMEnergyPar.hh"
#include "Offline/DbTables/inc/STMPedestals.hh"
#include "Offline/STMConditions/inc/STMEnergyCalib.hh"
#include "Offline/STMConfig/inc/STMEnergyCalibConfig.hh"

namespace mu2e {

class STMEnergyCalibMaker {
 public:
  STMEnergyCalibMaker(STMEnergyCalibConfig const& config) : _config(config) {}

  STMEnergyCalib::ptr_t fromFcl();
  STMEnergyCalib::ptr_t fromDb(STMEnergyPar::cptr_t sep_p,
                               STMPedestals::cptr_t ped_p);

 private:
  // this object needs to be thread safe,
  // _config should only be initialized once
  const STMEnergyCalibConfig _config;
};
}  // namespace mu2e

#endif
