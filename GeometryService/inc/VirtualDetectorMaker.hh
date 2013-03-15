#ifndef VirtualDetectorGeom_VirtualDetectorMaker_hh
#define VirtualDetectorGeom_VirtualDetectorMaker_hh
//
// Construct and return an VirtualDetector.
//
// $Id: VirtualDetectorMaker.hh,v 1.3 2013/03/15 15:52:04 kutschke Exp $
// $Author: kutschke $
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
