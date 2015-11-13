#ifndef GeometryService_DetectorSolenoidMaker_hh
#define GeometryService_DetectorSolenoidMaker_hh

#include <memory>

namespace mu2e  { class SimpleConfig; }
namespace mu2e  { class DetectorSolenoid; }
namespace mu2e  { class Beamline; }

namespace mu2e {
  class DetectorSolenoidMaker {
  public:
    static std::unique_ptr<DetectorSolenoid> make(const SimpleConfig& config, const Beamline& bl);
  };
}

#endif/* GeometryService_DetectorSolenoidMaker_hh */
