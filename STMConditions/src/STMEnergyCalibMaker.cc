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

  STMEnergyCalib::PedestalMap pmap;
  for (auto const& entry : _config.pedestals()) {
    std::string name = std::get<0>(entry);
    auto channel = STMChannel::findByName(name);
    if (!channel.isValid()) {
      throw cet::exception("STMENERGYCALIBMAKER_BAD_CHANNEL")
          << "STMEnergyCalibMaker::fromFcl called with bad channel name: "
          << name << "\n";
    }
    pmap[channel] = std::get<1>(entry);
  }

  auto ptr = std::make_shared<STMEnergyCalib>(cmap, pmap);
  return ptr;

}  // end fromFcl

STMEnergyCalib::ptr_t STMEnergyCalibMaker::fromDb(STMEnergyPar::cptr_t sep_p,
                                                  STMPedestals::cptr_t ped_p) {
  if (_config.verbose()) {
    std::cout << "STMEnergyCalibMaker::fromDb making STMEnergyCalib\n";
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

  STMEnergyCalib::PedestalMap pmap;
  for (auto const& row : ped_p->rows()) {
    pmap[row.channel()] = row.pedestal();
    if (_config.verbose()) {
      std::cout << std::setw(10) << row.channel().name() << std::fixed
                << std::setprecision(5) << std::setw(9) << row.pedestal()
                << "\n";
    }
  }

  auto ptr = std::make_shared<STMEnergyCalib>(cmap, pmap);
  return ptr;

}  // end fromDb

}  // namespace mu2e
