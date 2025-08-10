// clang-format off
// C++ includes
#include <iostream>

// Mu2e includes
#include "Offline/TrackerConditions/inc/TrkPanelMapEntity.hh"
#include <sstream>

using namespace std;

namespace mu2e {

  TrkPanelMapEntity::TrkPanelMapEntity(): ProditionsEntity(cxname) {
    _map.reserve(kMaxPanels);
  }
  
  void TrkPanelMapEntity::print(ostream& out) const {
  }

//-----------------------------------------------------------------------------
  void TrkPanelMapEntity::add(const TrkPanelMap::Row& Tpm) {
    _map.push_back(Tpm);
    const TrkPanelMap::Row* r = &_map.back();
    _tpm_by_mnid   [r->mnid()]              = r;
    _tpm_by_offline[r->plane()][r->panel()] = r;
    _tpm_by_online [r->dtc()  ][r->link() ] = r;
  }
} // namespace mu2e
