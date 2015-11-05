#ifndef GeometryService_TsdAMaker_hh
#define GeometryService_TsdAMaker_hh

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

#endif/* GeometryService_TsdAMaker_hh */
