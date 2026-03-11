// clang-format off
// C++ includes
#include <iostream>

// Mu2e includes
#include "Offline/TrackerConditions/inc/TrackerPanelMap.hh"
#include <sstream>

using namespace std;

namespace mu2e {

  TrackerPanelMap::TrackerPanelMap(): ProditionsEntity(cxname) {
    _map.reserve(kMaxPanels);
    for (int i=0; i<kMaxPanels; i++) {
      _tpm_by_mnid[i] = nullptr;
    }
    for (int plane=0; plane<kMaxPlanes; plane++) {
      for (int panel=0; panel<6; panel++) {
        _tpm_by_offline[plane][panel] = nullptr;
      }
    }
//-----------------------------------------------------------------------------
// N(links/DTC) = N(panels/plane) = 6 ... 1 DTC/plane
//-----------------------------------------------------------------------------
    for (int dtc=0; dtc<kMaxPlanes; dtc++) {
      for (int link=0; link<6; link++) {
        _tpm_by_online[dtc][link] = nullptr;
      }
    }
  }

  void TrackerPanelMap::print(ostream& out) const {
  }

//-----------------------------------------------------------------------------
  void TrackerPanelMap::add(const TrkPanelMap::Row& Tpm) {
    _map.push_back(Tpm);
    const TrkPanelMap::Row* r = &_map.back();
    _tpm_by_mnid   [r->mnid()       ]             = r;
    _tpm_by_offline[r->uniquePlane()][r->panel()] = r;
    _tpm_by_online [r->dtc()        ][r->link() ] = r;
  }
} // namespace mu2e
