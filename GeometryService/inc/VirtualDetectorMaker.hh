#ifndef VirtualDetectorGeom_VirtualDetectorMaker_hh
#define VirtualDetectorGeom_VirtualDetectorMaker_hh
//
// Construct and return an VirtualDetector.
//
// $Id: VirtualDetectorMaker.hh,v 1.2 2012/03/30 19:18:03 gandr Exp $
// $Author: gandr $
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
    static std::auto_ptr<VirtualDetector> make(const SimpleConfig& config);
  };

}  //namespace mu2e

#endif /* VirtualDetectorGeom_VirtualDetectorMaker_hh */
