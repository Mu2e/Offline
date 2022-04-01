#include "Offline/CaloConditions/inc/CalEnergyCalibMaker.hh"
#include "Offline/DataProducts/inc/CaloId.hh"
#include <vector>

/*
Written by S. Middleton
FIXME - placeholder
*/
namespace mu2e {
  typedef std::shared_ptr<CalEnergyCalib> ptr_t;

  ptr_t CalEnergyCalibMaker::fromFcl() {
    auto ptr = std::make_shared<CalEnergyCalib>(_config.roid());
    return ptr;

  } // end fromFcl

  ptr_t CalEnergyCalibMaker::fromDb(CalEnergyCalib::cptr_t) {
    // initially fill from fcl to get all the constants
    auto ptr = fromFcl();
    return ptr;
  } // end fromDb

}
