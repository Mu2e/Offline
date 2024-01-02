#include "Offline/CRVConditions/inc/CRVScintYieldMaker.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <iostream>

using namespace std;

namespace mu2e {

//***************************************************

CRVScintYield::ptr_t CRVScintYieldMaker::fromFcl() {
  if (_config.verbose()) {
    cout << "scintillation yield spread set to 0 for all CRV counters, because database is not used.\n";
  }

  size_t nCounters = GeomHandle<CosmicRayShield>()->getAllCRSScintillatorBars().size();

  CRVScintYield::ScintYieldVec svec(nCounters, 0.0);

  auto ptr = make_shared<CRVScintYield>(svec);
  return ptr;

}  // end fromFcl

//***************************************************

CRVScintYield::ptr_t CRVScintYieldMaker::fromDb(CRVScint::cptr_t sci_p) {
  if (_config.verbose()) {
    cout << "CRVScintYieldMaker::fromDb making CRVScintYield\n";
  }

  size_t nCounters = GeomHandle<CosmicRayShield>()->getAllCRSScintillatorBars().size();

  if (_config.verbose()) {
    cout << "CRVScintYieldMaker::fromDb checking for " << nCounters << " counters\n";
  }

  // require the db tables are the same length as geometry
  if (sci_p->nrow() != nCounters) {
    throw cet::exception("CRVSCINTYIELDMAKE_BAD_N_COUNTERS")
        << "CRVScintYieldMaker::fromDb bad counter counts: "
        << "  geometry: " << nCounters << "  CRVScint: " << sci_p->nrow() << "\n";
  }

  CRVScintYield::ScintYieldVec svec(nCounters, 0.0);

  for (auto const& row : sci_p->rows()) {
    svec[row.counter()] = row.scintYieldDeviation();
    if (_config.verbose()) {
      if (_config.verbose() > 1 || row.counter() < 5 ||
          CRVId::nBars - row.counter() <= 5) {
        cout << setw(10) << row.counter() << fixed << setprecision(3)
             << setw(10) << row.scintYieldDeviation() << "\n";
      }
    }
  }

  auto ptr = make_shared<CRVScintYield>(svec);
  return ptr;

}  // end fromDb

}  // namespace mu2e
