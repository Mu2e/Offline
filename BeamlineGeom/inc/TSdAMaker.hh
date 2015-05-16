#ifndef BEAMLINEGEOM_TSDAMAKER_HH
#define BEAMLINEGEOM_TSDAMAKER_HH

#include <memory>

namespace mu2e  { class SimpleConfig; }
namespace mu2e  { class TSdA; }
namespace mu2e  { class DetectorSolenoid; } 

namespace mu2e {
  class TSdAMaker {
  public:
    static std::unique_ptr<TSdA> make(const SimpleConfig& config, const DetectorSolenoid& ds);
  };
}

#endif/*BEAMLINEGEOM_TSDAMAKER_HH*/
