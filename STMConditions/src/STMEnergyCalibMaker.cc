#include "Offline/STMConditions/inc/STMEnergyCalibMaker.hh"
#include "cetlib_except/exception.h"
#include "TMath.h"
#include <iomanip>
#include <iostream>

namespace mu2e {

STMEnergyCalib::ptr_t STMEnergyCalibMaker::fromFcl() {
  if (_config.verbose()) {
    std::cout << "STMEnergyCalibMaker::fromFcl making nominal STMEnergyCalib\n";
  }

  STMEnergyCorr nominal{0.0, 1.0, 0.0};
  STMEnergyCalib::CalibMap cmap;

  cmap[STMChannel(STMChannel::enum_type::HPGe)] = nominal;
  cmap[STMChannel(STMChannel::enum_type::LaBr)] = nominal;

  auto ptr = std::make_shared<STMEnergyCalib>(cmap);
  return ptr;

}  // end fromFcl

STMEnergyCalib::ptr_t STMEnergyCalibMaker::fromDb(STMEnergyPar::cptr_t sep_p) {
  if (_config.verbose()) {
    std::cout << "STMEnergyCalibMaker::fromDb making nominal STMEnergyCalib\n";
  }

  STMEnergyCalib::CalibMap cmap;

  STMEnergyCorr corr;
  for (auto const& row : sep_p->rows()) {
    corr.p0 = row.p0();
    corr.p1 = row.p1();
    corr.p2 = row.p2();
    cmap[row.channel()] = corr;
    if (_config.verbose()) {
      std::cout << std::setw(10) << row.channel().name() << std::fixed
                << std::setprecision(5) << std::setw(9) << corr.p0
                << std::setw(9) << corr.p1 << std::setw(9) << corr.p2 << "\n";
    }
  }

  auto ptr = std::make_shared<STMEnergyCalib>(cmap);
  return ptr;

}  // end fromDb

}  // namespace mu2e
