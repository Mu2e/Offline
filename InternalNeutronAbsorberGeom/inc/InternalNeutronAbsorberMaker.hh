#ifndef INTERNALNEUTRONABSORBERMAKER_HH
#define INTERNALNEUTRONABSORBERMAKER_HH

#include <memory>

namespace mu2e  { class SimpleConfig; }
namespace mu2e  { class InternalNeutronAbsorber; }
namespace mu2e  { class DetectorSolenoid; } 

namespace mu2e {
  class InternalNeutronAbsorberMaker {
  public:
    static std::unique_ptr<InternalNeutronAbsorber> make(const SimpleConfig& config, const DetectorSolenoid& ds);
  };
}

#endif/*INTERNALNEUTRONABSORBERMAKER_HH*/
