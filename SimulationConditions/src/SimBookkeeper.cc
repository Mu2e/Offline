// Mu2e includes
#include "Offline/SimulationConditions/inc/SimBookkeeper.hh"
#include "cetlib_except/exception.h"

namespace mu2e {

  double SimBookkeeper::getEff(const std::string& name) const {
    auto it = _effs.find(name);
    if (it == _effs.end()) {
      cet::exception ex("SIMBOOKKEEPER_MISSING_TAG");
      ex << "SimBookkeeper has no efficiency with tag \"" << name
         << "\". Available tags:";
      for (const auto& i_eff : _effs) {
        ex << " " << i_eff.first;
      }
      ex << "\n";
      throw ex;
    }
    return it->second;
  }

  const std::string SimBookkeeper::print() const {
    std::stringstream out;
    print(out);
    return out.str();
  }

  void SimBookkeeper::print(std::ostream& os) const {
    os << "Efficiencies in " << name() << ":" << std::endl;
    for (const auto& i_eff : _effs) {
      os << i_eff.first << " = " << i_eff.second << std::endl;
    }
    os << std::endl;
  }
}
