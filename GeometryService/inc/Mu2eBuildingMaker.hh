#ifndef MU2EBUILDINGMAKER_HH
#define MU2EBUILDINGMAKER_HH

#include <memory>

namespace mu2e { class SimpleConfig; }
namespace mu2e { class Mu2eBuilding; }

namespace mu2e {
  class Mu2eBuildingMaker {
    std::auto_ptr<Mu2eBuilding> _b;
  public:
    explicit Mu2eBuildingMaker(const SimpleConfig& config);
    
    // interface to GeometryService
    std::auto_ptr<Mu2eBuilding> getPtr() { return _b; }
  };
}

#endif/*MU2EBUILDINGMAKER_HH*/
