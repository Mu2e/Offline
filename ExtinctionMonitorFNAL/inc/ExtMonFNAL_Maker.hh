#ifndef EXTMONFNAL_MAKER_HH
#define EXTMONFNAL_MAKER_HH

#include <memory>

namespace mu2e { class SimpleConfig; }
namespace mu2e { namespace ExtMonFNAL { class ExtMon; } }

namespace mu2e {
  namespace ExtMonFNAL {

    class ExtMonMaker {
      std::auto_ptr<ExtMon> m_det;
    public:
      ExtMonMaker(const SimpleConfig& config);

      // interface to GeometryService
      std::auto_ptr<ExtMon> getDetectorPtr() { return m_det; }
    };
  }
}

#endif/*EXTMONFNAL_MAKER_HH*/
