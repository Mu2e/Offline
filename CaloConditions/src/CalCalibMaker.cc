#include "Offline/CaloConditions/inc/CalCalibMaker.hh"
#include <vector>
using namespace std;
/*
Written by S. Middleton - imports calibration information for the calo from either fcl or the database.
*/
namespace mu2e {
  typedef std::shared_ptr<CalCalib> ptr_t;

/*******************************************************************/
  ptr_t CalCalibMaker::fromFcl() {

  if (_config.verbose()) {
    std::cout << "CalCalibMaker::fromFcl making nominal CalCalib\n";
  }
  CalCalibPar nominal(_config.ADC2MeV(), _config.timeoffset());

  size_t nChan = CaloConst::_nChannel;

  if (_config.verbose()) {
    std::cout << "CalCalibMaker::fromFcl filling " << nChan << " channels\n";
    std::cout << "CalCalibMaker::fromFcl nominal " << fixed << setprecision(3)
         << setw(10) << nominal.ADC2MeV() << setprecision(3) << setw(10)
         << nominal.timeOffset() << setprecision(3) << setw(10) << "\n";
  }

  CalCalib::CalibVec cvec(nChan, nominal);
  auto ptr = std::make_shared<CalCalib>(cvec);
  return ptr;

  } // end fromFcl
/*******************************************************************/


/*******************************************************************/

  ptr_t CalCalibMaker::fromDb(const CalEnergyCalib& ecalib,
                              const CalTimeCalib& tcalib) {
    // initially fill from fcl to get all the constants
    if (_config.verbose()) {
      cout << "CalCalibMaker::fromDb making CalCalib\n";
    }

    size_t nChan = CaloConst::_nChannel;
    size_t nChanDB = CaloConst::_nChannelDB;

    if (_config.verbose()) {
      cout << "CalCalibMaker::fromDb checking for " << nChanDB << " channels\n";
    }

    // check the db tables length (TODO remove?)
    if (ecalib.nrow() != nChanDB || tcalib.nrow() != nChanDB) {
      throw cet::exception("CALCALIBMAKE_BAD_N_CHANNEL")
      << "CalCalibMaker::fromDb bad channel counts: "
      << "  geometry: " << nChan
      << "  Nenergy: " << ecalib.nrow()
      << "  Ntime: " << ecalib.nrow()
      << "\n";
    }

    CalCalib::CalibVec cvec;

    //Loop up to _nChannel (skip the spare DB channels)
    for (CaloConst::CaloSiPMId_type ind=0; ind<nChan; ind++) {
      auto roid = CaloSiPMId(ind);
      cvec.emplace_back(ecalib.row(roid).ADC2MeV(), tcalib.row(roid).tcorr());
    }

    auto ptr = make_shared<CalCalib>(cvec);
    return ptr;

  } // end fromDb
/*******************************************************************/
}
