#ifndef ValTrkCrvSectorAxis_HH_
#define ValTrkCrvSectorAxis_HH_

#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "TAxis.h"
#include <map>
#include <string>
#include <vector>

namespace mu2e {

// fixed, labelled histogram axis of the CRV sector SurfaceIds, so the per-sector histograms do not depend on the job
// geometry; the last bin collects sectors without a SurfaceId
class ValTrkCrvSectorAxis {
 public:
  ValTrkCrvSectorAxis() {
    for (auto const& [id, name] : SurfaceIdDetail::names()) {
      if (name.rfind("CRV_", 0) != 0 || id == SurfaceIdDetail::CRV_StrongBack) continue;
      _byName[name] = _bin[id] = static_cast<int>(_names.size());
      _names.push_back(name);
    }
    _names.push_back("unknown");
  }
  int n() const { return static_cast<int>(_names.size()); }
  int unknown() const { return n() - 1; }
  double lo() const { return -0.5; }
  double hi() const { return n() - 0.5; }
  int bin(SurfaceId const& sid) const {
    auto it = _bin.find(sid.id().id());
    return it == _bin.end() ? unknown() : it->second;
  }
  int bin(std::string const& name) const {
    auto it = _byName.find(name);
    return it == _byName.end() ? unknown() : it->second;
  }
  void label(TAxis* ax) const {
    for (int i = 0; i < n(); ++i) ax->SetBinLabel(i + 1, _names[i].c_str());
  }

 private:
  std::map<SurfaceIdDetail::enum_type, int> _bin;
  std::map<std::string, int> _byName;
  std::vector<std::string> _names;
};
}  // namespace mu2e

#endif
