#include "Offline/CRVConditions/inc/CRVStatusMaker.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <iostream>

using namespace std;

namespace mu2e {

//***************************************************

CRVStatus::ptr_t CRVStatusMaker::fromFcl() {
  if (_config.verbose()) {
    cout << "CRVStatusMaker::fromFcl making nominal CRVStatus\n";
  }

  // this is a fcl-derived vector of tuples (channel,status)
  auto fvec = _config.status();

  CRVStatus::StatusMap smap;

  for (auto const& ss : fvec) {
    uint16_t channel = get<0>(ss);
    int stat = get<1>(ss);
    smap[channel] = stat;
    if (_config.verbose()) {
      cout << "   CRVStatusMaker loading fcl  channel " << setw(6) << channel
           << " with status " << setw(6) << stat << "\n";
    }
  }

  auto ptr = make_shared<CRVStatus>(smap);
  return ptr;

}  // end fromFcl

//***************************************************

CRVStatus::ptr_t CRVStatusMaker::fromDb(CRVBadChan::cptr_t cbc_p) {
  if (_config.verbose()) {
    cout << "CRVStatusMaker::fromDb making CRVStatus\n";
  }

  CRVStatus::StatusMap smap;

  for (auto const& row : cbc_p->rows()) {
    smap[row.channel()] = row.status();
    if (_config.verbose()) {
      cout << "   CRVStatusMaker loading channel " << setw(6) << row.channel()
           << " with status " << setw(6) << row.status() << "\n";
    }
  }

  auto ptr = make_shared<CRVStatus>(smap);
  return ptr;

}  // end fromDb

}  // namespace mu2e
