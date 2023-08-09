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
  CalCalibPar nominal(_config.ADC2MeV(), _config.ECombAlgID(), _config.timeoffset());

  size_t nChan = CaloConst::_nChannel;

  if (_config.verbose()) {
    std::cout << "CalCalibMaker::fromFcl filling " << nChan << " channels\n";
    std::cout << "CalCalibMaker::fromFcl nominal " << fixed << setprecision(3)
         << setw(10) << nominal.ADC2MeV() << setprecision(3) << setw(10)
         << nominal.ECombAlgID() << setprecision(3) << setw(10)
         << nominal.timeOffset() << setprecision(3) << setw(10) << "\n";
  }

  CalCalib::CalibVec cvec(nChan, nominal);
  auto ptr = std::make_shared<CalCalib>(cvec);
  return ptr;

  } // end fromFcl
/*******************************************************************/


/*******************************************************************/

  ptr_t CalCalibMaker::fromDb(CalEnergyCalib::cptr_t ecalib) { //TODO - input from time table too
    // initially fill from fcl to get all the constants
    if (_config.verbose()) {
      cout << "CalCalibMaker::fromDb making CalCalib\n";
    }

    size_t nChan = CaloConst::_nChannel;

    if (_config.verbose()) {
      cout << "CalCalibMaker::fromDb checking for " << nChan << " channels\n";
    }

  // require the db tables are the same length as geometry
  if (ecalib->nrow() != nChan ) {
    throw cet::exception("CALCALIBMAKE_BAD_N_CHANNEL")
        << "CalCalibMaker::fromDb bad channel counts: "
        << "  geometry: " << nChan << "  CalSiPM: " << ecalib->nrow()<< "\n";
  }

  CalCalib::CalibVec cvec(nChan, CalCalibPar(0.0, 0.0, 0.0));

  for (auto const& row : ecalib->rows()) {
    cvec[row.roid().id()] = CalCalibPar(row.ADC2MeV(), row.algID(),0); //TODO - time offset needs setting from time table
  }

  auto ptr = make_shared<CalCalib>(cvec);
  return ptr;
    //auto ptr = fromFcl();
    //return ptr;

  } // end fromDb
/*******************************************************************/
}
