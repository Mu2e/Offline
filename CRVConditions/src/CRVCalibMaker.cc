#include "Offline/CRVConditions/inc/CRVCalibMaker.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <iostream>

using namespace std;

namespace mu2e {

//***************************************************

CRVCalib::ptr_t CRVCalibMaker::fromFcl() {
  if (_config.verbose()) {
    cout << "CRVCalibMaker::fromFcl making nominal CRVCalib\n";
  }

  CRVCalibPar nominal(_config.pedestal(), _config.pulseHeight(),
                      _config.pulseArea(), _config.timeOffset());

  size_t nChan =
      GeomHandle<CosmicRayShield>()->getAllCRSScintillatorBars().size() *
      CRVId::nChanPerBar;

  if (_config.verbose()) {
    cout << "CRVCalibMaker::fromFcl filling " << nChan << " channels\n";
    cout << "CRVCalibMaker::fromFcl nominal " << fixed << setprecision(3)
         << setw(8) << nominal.pedestal() << setprecision(3) << setw(8)
         << nominal.pulseHeight() << setprecision(3) << setw(8)
         << nominal.pulseArea() << setprecision(3) << setw(8)
         << nominal.timeOffset() << "\n";
  }

  CRVCalib::CalibVec cvec(nChan, nominal);

  auto ptr = make_shared<CRVCalib>(cvec);
  return ptr;

}  // end fromFcl

//***************************************************

CRVCalib::ptr_t CRVCalibMaker::fromDb(CRVSiPM::cptr_t sip_p,
                                      CRVTime::cptr_t tim_p) {
  if (_config.verbose()) {
    cout << "CRVCalibMaker::fromDb making CRVCalib\n";
  }

  size_t nChan =
      GeomHandle<CosmicRayShield>()->getAllCRSScintillatorBars().size() *
      CRVId::nChanPerBar;

  if (_config.verbose()) {
    cout << "CRVCalibMaker::fromDb expecting " << nChan << " channels\n";
  }

  // require the db tables are the same length as geometry
  if (sip_p->nrow() != nChan || tim_p->nrow() != nChan) {
    throw cet::exception("CRVCALIBMAKE_BAD_N_CHANNEL")
        << "CRVCalibMaker::fromDb bad channel counts: "
        << "  geometry: " << nChan << "  CRVSiPM: " << sip_p->nrow()
        << "  CRVTime: " << tim_p->nrow() << "\n";
  }

  CRVCalib::CalibVec cvec(nChan, CRVCalibPar(0.0, 0.0, 0.0, 0.0));

  for (auto const& row : sip_p->rows()) {
    float timeOffset = tim_p->row(row.channel()).timeOffset();
    cvec[row.channel()] = CRVCalibPar(row.pedestal(), row.pulseHeight(),
                                      row.pulseArea(), timeOffset);
    if (_config.verbose()) {
      if (_config.verbose() > 1 || row.channel() < 5 ||
          CRVId::nChannels - row.channel() <= 5) {
        cout << setw(10) << row.channel() << fixed << setprecision(3) << setw(8)
             << row.pedestal() << setprecision(3) << setw(8)
             << row.pulseHeight() << setprecision(3) << setw(8)
             << row.pulseArea() << setprecision(3) << setw(8) << timeOffset
             << "\n";
      }
    }
  }

  auto ptr = make_shared<CRVCalib>(cvec);
  return ptr;

}  // end fromDb

}  // namespace mu2e
