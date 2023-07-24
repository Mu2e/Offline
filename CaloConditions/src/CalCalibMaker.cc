#include "Offline/CaloConditions/inc/CalCalibMaker.hh"
#include <vector>

/*
Written by S. Middleton
*/
namespace mu2e {
  typedef std::shared_ptr<CalCalib> ptr_t;

/*******************************************************************/
  ptr_t CalCalibMaker::fromFcl() {
  
  if (_config.verbose()) {
    cout << "CRalCalibMaker::fromFcl making nominal CalCalib\n";
  }
  CalCalibPar nominal(_config.ADC2MeV(), _config.timeoffset(), _config.algID());

  size_t nChan = CaloConst::_nChannel;

  if (_config.verbose()) {
    cout << "CalCalibMaker::fromFcl filling " << nChan << " channels\n";
    cout << "CalCalibMaker::fromFcl nominal " << fixed << setprecision(3)
         << setw(10) << nominal.ADC2MeV() << setprecision(3) << setw(10)
         << nominal.timeoffset() << setprecision(3) << setw(10)
         << nominal.algID() << setprecision(3) << setw(10) << "\n";
  }
  
  CalCalib::CalibVec cvec(nChan, nominal);
  auto ptr = std::make_shared<CalCalib>(cvec);
  return ptr;

  } // end fromFcl
/*******************************************************************/


/*******************************************************************/

  ptr_t CalCalibMaker::fromDb(CalEnergyCalib::cptr_t) {
    // initially fill from fcl to get all the constants

    auto ptr = fromFcl();
    return ptr;

  } // end fromDb
/*******************************************************************/
}
