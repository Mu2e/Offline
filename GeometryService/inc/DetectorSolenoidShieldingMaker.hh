#ifndef DETECTORSOLENOIDSHIELDINGMAKER_HH
#define DETECTORSOLENOIDSHIELDINGMAKER_HH

#include <memory>

namespace mu2e  { class SimpleConfig;              }
namespace mu2e  { class DetectorSolenoidShielding; }
namespace mu2e  { class DetectorSolenoid;          }

namespace mu2e {
  class DetectorSolenoidShieldingMaker {
  public:
    static std::unique_ptr<DetectorSolenoidShielding> make(const SimpleConfig& config, const DetectorSolenoid& ds);
  };
}

#endif/*DETECTORSOLENOIDSHIELDINGMAKER_HH*/
