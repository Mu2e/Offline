#ifndef EXTMONFNAL_MAKER_HH
#define EXTMONFNAL_MAKER_HH

#include <memory>
#include <string>

#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALSensorStack.hh"

namespace mu2e { class SimpleConfig; }
namespace mu2e { namespace ExtMonFNAL { class ExtMon; } }
namespace mu2e { class ExtMonFNALBuilding; }
namespace CLHEP { class Hep3Vector; }
namespace CLHEP { class Hep3Rotation; }


namespace mu2e {
  namespace ExtMonFNAL {

    class ExtMonMaker {
      static ExtMonFNALSensorStack readStack(const SimpleConfig& c,
                                             const std::string& prefix,
                                             const CLHEP::Hep3Vector& refPointInMu2e,
                                             const CLHEP::HepRotation& rotationInMu2e
                                             );

    public:
      static std::unique_ptr<ExtMon> make(const SimpleConfig& config, const ExtMonFNALBuilding& room);
    };
  }
}

#endif/*EXTMONFNAL_MAKER_HH*/
