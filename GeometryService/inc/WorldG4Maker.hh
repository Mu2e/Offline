#ifndef WORLDG4MAKER_HH
#define WORLDG4MAKER_HH

#include <memory>

namespace mu2e { class SimpleConfig; }
namespace mu2e { class WorldG4; }

namespace mu2e {
  class WorldG4Maker {
    std::auto_ptr<WorldG4> _wg4;
  public:
    explicit WorldG4Maker(const SimpleConfig& config);
    
    // interface to GeometryService
    std::auto_ptr<WorldG4> getPtr() { return _wg4; }
  };
}

#endif/*WORLDG4MAKER_HH*/
