#ifndef DETECTORSOLENOIDMAKER_HH
#define DETECTORSOLENOIDMAKER_HH

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

#endif/*DETECTORSOLENOIDMAKER_HH*/
