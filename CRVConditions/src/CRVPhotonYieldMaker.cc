#include "Offline/CRVConditions/inc/CRVPhotonYieldMaker.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "cetlib_except/exception.h"
#include <iomanip>
#include <iostream>

using namespace std;

namespace mu2e {

//***************************************************

CRVPhotonYield::ptr_t CRVPhotonYieldMaker::fromFcl() {
  if (_config.verbose()) {
    cout << "photon yield deviations set to 0 for all CRV channels, because database is not used.\n";
  }

  size_t nChan =
      GeomHandle<CosmicRayShield>()->getAllCRSScintillatorBars().size() *
      CRVId::nChanPerBar;

  CRVPhotonYield::PhotonYieldVec svec(nChan, 0.0);

  auto ptr = make_shared<CRVPhotonYield>(svec);
  return ptr;

}  // end fromFcl

//***************************************************

CRVPhotonYield::ptr_t CRVPhotonYieldMaker::fromDb(CRVPhoton::cptr_t sci_p) {
  if (_config.verbose()) {
    cout << "CRVPhotonYieldMaker::fromDb making CRVPhotonYield\n";
  }

  size_t nChan =
      GeomHandle<CosmicRayShield>()->getAllCRSScintillatorBars().size() *
      CRVId::nChanPerBar;

  if (_config.verbose()) {
    cout << "CRVPhotonYieldMaker::fromDb checking for " << nChan << " channels\n";
  }

  // require the db tables are the same length as geometry
  if (sci_p->nrow() != nChan) {
    throw cet::exception("CRVPHOTONYIELDMAKE_BAD_N_COUNTERS")
        << "CRVPhotonYieldMaker::fromDb bad channel counts: "
        << "  geometry: " << nChan << "  CRVPhoton: " << sci_p->nrow() << "\n";
  }

  CRVPhotonYield::PhotonYieldVec svec(nChan, 0.0);

  for (auto const& row : sci_p->rows()) {
    svec[row.channel()] = row.photonYieldDeviation();
    if (_config.verbose()) {
      if (_config.verbose() > 1 || row.channel() < 5 ||
          CRVId::nChannels - row.channel() <= 5) {
        cout << setw(10) << row.channel() << fixed << setprecision(3)
             << setw(10) << row.photonYieldDeviation() << "\n";
      }
    }
  }

  auto ptr = make_shared<CRVPhotonYield>(svec);
  return ptr;

}  // end fromDb

}  // namespace mu2e
