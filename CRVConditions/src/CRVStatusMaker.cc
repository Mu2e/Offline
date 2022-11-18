#include "Offline/CRVConditions/inc/CRVStatusMaker.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <iostream>

namespace mu2e {

CRVStatus::ptr_t CRVStatusMaker::fromFcl() {
  if (_config.verbose()) {
    std::cout << "CRVStatusMaker::fromFcl making nominal CRVStatus\n";
  }

  // this is a fcl vector of tuples (channel,status)
  auto fvec = _config.status();

  CRVStatus::StatusMap smap;

  for (auto const& ss : fvec) {
    std::size_t channel = std::get<0>(ss);
    int stat = std::get<1>(ss);
    smap[channel] = stat;
    if (_config.verbose()) {
      std::cout << "   CRVStatusMaker loading   channel " << std::setw(6)
                << channel << " with status " << std::fixed
                << std::setprecision(3) << std::setw(9) << stat << "\n";
    }
  }

  auto ptr = std::make_shared<CRVStatus>(smap);
  return ptr;

}  // end fromFcl

CRVStatus::ptr_t CRVStatusMaker::fromDb(CRVBadChan::cptr_t cbc_p) {
  if (_config.verbose()) {
    std::cout << "CRVStatusMaker::fromDb making CRVStatus\n";
  }

  CRVStatus::StatusMap smap;

  for (auto const& row : cbc_p->rows()) {
    smap[row.channel()] = row.status();
    if (_config.verbose()) {
      std::cout << "   CRVStatusMaker loading   channel " << std::setw(6)
                << row.channel() << " with status " << std::fixed
                << std::setprecision(3) << std::setw(9) << row.status() << "\n";
    }
  }

  auto ptr = std::make_shared<CRVStatus>(smap);
  return ptr;

}  // end fromDb

}  // namespace mu2e
