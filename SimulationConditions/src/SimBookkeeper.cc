// Mu2e includes
#include "SimulationConditions/inc/SimBookkeeper.hh"

namespace mu2e {

  const std::string SimBookkeeper::print() const {
    std::stringstream out;
    print(out);
    return out.str();
  }

  void SimBookkeeper::print(std::ostream& os) const {
    os << "Efficiencies in " << _name << ":" << std::endl;
    for (const auto& i_eff : _effs) {
      os << i_eff.first << " = " << i_eff.second << std::endl;
    }
    os << std::endl;
  }
}
