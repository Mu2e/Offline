#ifndef ProductionSolenoidGeom_PSEnclosureMaker_hh
#define ProductionSolenoidGeom_PSEnclosureMaker_hh
//
// Class to construct and return PSEnclosure
//
//
// Original author Andrei Gaponenko
//

#include <string>
#include <memory>

namespace CLHEP { class Hep3Vector; }

namespace mu2e {

  class PSEnclosure;
  class SimpleConfig;

  class PSEnclosureMaker {
  public:

    static std::unique_ptr<PSEnclosure>  make(const SimpleConfig& config,

                                            // The center of the downstream surface of the PS
                                            const CLHEP::Hep3Vector& psEndRefPoint);
  private:
    static void readWindow(int i, const SimpleConfig& c, const CLHEP::Hep3Vector& winRefPoint, std::unique_ptr<PSEnclosure>& res);

  };
}  //namespace mu2e

#endif /* ProductionSolenoidGeom_PSEnclosureMaker_hh */
