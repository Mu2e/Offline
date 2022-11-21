#include "Offline/CRVConditions/inc/CRVCalibMaker.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <iostream>

namespace mu2e {

CRVCalib::ptr_t CRVCalibMaker::fromFcl() {
  if (_config.verbose()) {
    std::cout << "CRVCalibMaker::fromFcl making nominal CRVCalib\n";
  }

  CRVCalibPar nominal(_config.pedestal(), _config.height(), _config.area(),
                      _config.timeOffset());
  CRVCalib::CalibVec cvec(CRVId::nChannels, nominal);

  auto ptr = std::make_shared<CRVCalib>(cvec);
  return ptr;

}  // end fromFcl

CRVCalib::ptr_t CRVCalibMaker::fromDb(CRVSiPM::cptr_t sip_p,
                                      CRVTime::cptr_t tim_p) {
  if (_config.verbose()) {
    std::cout << "CRVCalibMaker::fromDb making CRVCalib\n";
  }

  CRVCalib::CalibVec cvec(CRVId::nChannels, {0.0, 0.0, 0.0, 0.0});

  for (auto const& row : sip_p->rows()) {
    float timeOffset = tim_p->row(row.channel()).timeOffset();
    cvec[row.channel()] =
        CRVCalibPar(row.pedestal(), row.height(), row.area(), timeOffset);
    if (_config.verbose()) {
      if (_config.verbose() > 1 || row.channel() < 5 ||
          CRVId::nChannels - row.channel() <= 5) {
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
