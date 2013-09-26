#ifndef MU2EBUILDINGMAKER_HH
#define MU2EBUILDINGMAKER_HH

#include <memory>

namespace mu2e {

  class SimpleConfig;
  class Mu2eBuilding;
  class BuildingBasics;
  class Beamline;
  class ProtonBeamDump;

  class Mu2eBuildingMaker {
  public:
    static std::unique_ptr<Mu2eBuilding> make(const SimpleConfig& config,
                                              const BuildingBasics& basics,
                                              const Beamline& bl,
                                              const ProtonBeamDump& dump);
  };
}

#endif/*MU2EBUILDINGMAKER_HH*/
