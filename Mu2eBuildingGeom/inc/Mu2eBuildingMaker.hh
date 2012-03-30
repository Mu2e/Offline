#ifndef MU2EBUILDINGMAKER_HH
#define MU2EBUILDINGMAKER_HH

#include <memory>

namespace mu2e { class SimpleConfig; }
namespace mu2e { class Mu2eBuilding; }
namespace mu2e { class ProtonBeamDump; }

namespace mu2e {
  class Mu2eBuildingMaker {
  public:
    static std::auto_ptr<Mu2eBuilding> make(const SimpleConfig& config, const ProtonBeamDump& dump);
  };
}

#endif/*MU2EBUILDINGMAKER_HH*/
