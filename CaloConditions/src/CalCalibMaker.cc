#include "Offline/CaloConditions/inc/CalCalibMaker.hh"
#include <vector>

/*
Written by S. Middleton
*/
namespace mu2e {
  typedef std::shared_ptr<CalCalib> ptr_t;

  ptr_t CalCalibMaker::fromFcl() {
  CalCalibPar nominal(_config.ADC2MeV(), _config.timeoffset(), _config.algID());

  size_t nChan = CaloConst::_nChannel;

  CalCalib::CalibVec cvec(nChan, nominal);
  auto ptr = std::make_shared<CalCalib>(cvec);
  return ptr;

  } // end fromFcl

  ptr_t CalCalibMaker::fromDb(CalEnergyCalib::cptr_t) {
    // initially fill from fcl to get all the constants

    auto ptr = fromFcl();
    return ptr;

  } // end fromDb

}
