#ifndef EXTMONFNAL_MAKER_HH
#define EXTMONFNAL_MAKER_HH

#include <memory>

namespace mu2e { class SimpleConfig; }
namespace mu2e { namespace ExtMonFNAL { class ExtMon; } }
namespace mu2e { class ProtonBeamDump; }

namespace mu2e {
  namespace ExtMonFNAL {

    class ExtMonMaker {
    public:
      static std::auto_ptr<ExtMon> make(const SimpleConfig& config, const ProtonBeamDump& dump);
    };
  }
}

#endif/*EXTMONFNAL_MAKER_HH*/
