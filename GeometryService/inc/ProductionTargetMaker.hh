#ifndef PRODUCTIONTARGETMAKER_HH
#define PRODUCTIONTARGETMAKER_HH

#include <memory>

namespace mu2e { class SimpleConfig; }
namespace mu2e { class ProductionTarget; }

namespace mu2e {
  class ProductionTargetMaker {
  public:
    static std::unique_ptr<ProductionTarget> make(const SimpleConfig& config, double solenoidOffset);
  };
}

#endif/*PRODUCTIONTARGETMAKER_HH*/
