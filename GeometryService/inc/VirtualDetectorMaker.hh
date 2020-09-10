#ifndef VirtualDetectorGeom_VirtualDetectorMaker_hh
#define VirtualDetectorGeom_VirtualDetectorMaker_hh
//
// Construct and return an VirtualDetector.
//
//
// Original author Peter Shanahan
//

#include <memory>
#include <string>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class SimpleConfig;
  class VirtualDetector;

  class VirtualDetectorMaker {
  public:
    static std::unique_ptr<VirtualDetector> make(const SimpleConfig& config);
  };

}  //namespace mu2e

#endif /* VirtualDetectorGeom_VirtualDetectorMaker_hh */
