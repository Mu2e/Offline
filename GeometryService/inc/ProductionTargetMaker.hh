#ifndef PRODUCTIONTARGETMAKER_HH
#define PRODUCTIONTARGETMAKER_HH

#include <memory>

namespace mu2e { class SimpleConfig; }
namespace mu2e { class ProductionTarget; }

namespace mu2e {
  class ProductionTargetMaker {
    std::auto_ptr<ProductionTarget> m_det;
  public:
    ProductionTargetMaker(const SimpleConfig& config);
    
    // interface to GeometryService
    std::auto_ptr<ProductionTarget> getDetectorPtr() { return m_det; }
  };
}

#endif/*PRODUCTIONTARGETMAKER_HH*/
