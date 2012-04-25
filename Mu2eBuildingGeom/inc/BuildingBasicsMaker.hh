#ifndef BUILDINGBASICSMAKER_HH
#define BUILDINGBASICSMAKER_HH

#include <memory>

namespace mu2e { class SimpleConfig; }
namespace mu2e { class BuildingBasics; }

namespace mu2e {
  class BuildingBasicsMaker {
  public:
    static std::auto_ptr<BuildingBasics> make(const SimpleConfig& config);
  };
}

#endif/*BUILDINGBASICSMAKER_HH*/
