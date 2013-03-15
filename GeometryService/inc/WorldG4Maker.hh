#ifndef WORLDG4MAKER_HH
#define WORLDG4MAKER_HH

#include <memory>

namespace mu2e { class SimpleConfig; }
namespace mu2e { class WorldG4; }

namespace mu2e {
  class WorldG4Maker {
  public:
    static std::unique_ptr<WorldG4> make(const SimpleConfig& config);
  };
}

#endif/*WORLDG4MAKER_HH*/
