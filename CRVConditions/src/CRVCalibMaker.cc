#include "Offline/CRVConditions/inc/CRVCalibMaker.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <iostream>

namespace mu2e {

std::size_t CRVCalibMaker::nChannels() const {
  // CRV has several geometries, so get size from there
  GeomHandle<CosmicRayShield> CRS;
  size_t nChan = CRS->getAllCRSScintillatorBars().size() * CRVId::nChanPerBar;
  return nChan;
}

CRVCalib::ptr_t CRVCalibMaker::fromFcl() {
  if (_config.verbose()) {
    std::cout << "CRVCalibMaker::fromFcl making nominal CRVCalib\n";
  }

  CRVCalibPar nominal(_config.pedestal(), _config.height(), _config.area(),
                      _config.timeOffset());

  // CRV has several geometries, so get size from there
  std::size_t nChan = nChannels();

  CRVCalib::CalibVec cvec(nChan, nominal);

  auto ptr = std::make_shared<CRVCalib>(cvec);
  return ptr;

}  // end fromFcl

CRVCalib::ptr_t CRVCalibMaker::fromDb(CRVSiPM::cptr_t sip_p,
                                      CRVTime::cptr_t tim_p) {
  if (_config.verbose()) {
    std::cout << "CRVCalibMaker::fromDb making CRVCalib\n";
  }

  // CRV has several geometries, so get size from there
  std::size_t nChan = nChannels();

  // require the db tables are the same length as geometry
  if (sip_p->nrow() != nChan || tim_p->nrow() != nChan) {
    throw cet::exception("CRVCALIBMAKE_BAD_N_CHANNEL")
        << "CRVCalibMaker::fromDb bad channel counts: "
        << "  geometry: " << nChan << "  CRVSiPM: " << sip_p->nrow()
        << "  CRVTime: " << tim_p->nrow() << "\n";
  }

  CRVCalib::CalibVec cvec(nChan, {0.0, 0.0, 0.0, 0.0});

  for (auto const& row : sip_p->rows()) {
    float timeOffset = tim_p->row(row.channel()).timeOffset();
    cvec[row.channel()] =
        CRVCalibPar(row.pedestal(), row.height(), row.area(), timeOffset);
    if (_config.verbose()) {
      if (_config.verbose() > 1 || row.channel() < 5 ||
          nChan - row.channel() <= 5) {
        std::cout << std::setw(10) << row.channel() << std::fixed
                  << std::setprecision(3) << std::setw(8) << row.pedestal()
                  << std::setprecision(3) << std::setw(8) << row.height()
                  << std::setprecision(3) << std::setw(8) << row.area()
                  << std::setprecision(3) << std::setw(8) << timeOffset << "\n";
      }
    }
  }

  auto ptr = std::make_shared<CRVCalib>(cvec);
  return ptr;

}  // end fromDb

}  // namespace mu2e
