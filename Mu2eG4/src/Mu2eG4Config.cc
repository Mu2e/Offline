#include "Mu2eG4/inc/Mu2eG4Config.hh"
#include "Mu2eG4/inc/Mu2eG4PrimaryType.hh"

namespace mu2e {
  namespace Mu2eG4Config {
    std::string Inputs_::primaryType_docstring() {
      std::string res{"The type of inputs used to create G4 primaries.  The allowed values are:"};
      for(const auto& s: Mu2eG4PrimaryType::all_names()) {
        res += " " + s;
      }
      return res;
    }
  }
}
