#include "Offline/CaloConditions/inc/CalEnergyCalibMaker.hh"
#include "Offline/DataProducts/inc/CaloId.hh"
#include <vector>

/*
Written by S. Middleton
FIXME - placeholder
*/
namespace mu2e {
  typedef std::shared_ptr<CalCalibConstant> ptr_t;

  ptr_t CalEnergyCalibMaker::fromFcl() {
    auto ptr = std::make_shared<CalCalibConstant>(_config.roid());
    return ptr;

  } // end fromFcl

  ptr_t CalEnergyCalibMaker::fromDb(CalCalibConstant::cptr_t) {
    // initially fill from fcl to get all the constants
    auto ptr = fromFcl();
    return ptr;
  } // end fromDb

}
