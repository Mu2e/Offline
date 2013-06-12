#ifndef EXTERNALNEUTRONABSORBERMAKER_HH
#define EXTERNALNEUTRONABSORBERMAKER_HH

#include <memory>

namespace mu2e  { class SimpleConfig; }
namespace mu2e  { class ExternalNeutronAbsorber; }
namespace mu2e  { class DetectorSolenoid; } 

namespace mu2e {
  class ExternalNeutronAbsorberMaker {
  public:
    static std::unique_ptr<ExternalNeutronAbsorber> make(const SimpleConfig& config, const DetectorSolenoid& ds);
  };
}

#endif/*EXTERNALNEUTRONABSORBERMAKER_HH*/
