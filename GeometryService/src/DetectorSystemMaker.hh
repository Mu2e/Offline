#ifndef GeometryService_src_DetectorSystemMaker_hh
#define GeometryService_src_DetectorSystemMaker_hh
//
// Construct a DetectorSystem object.
//
//
// Original author Rob Kutschke
//

// C++ includes
#include <memory>

// Mu2e includes.
#include "Offline/Mu2eInterfaces/inc/Detector.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"

namespace mu2e {

  // Forward references
  class SimpleConfig;

  class DetectorSystemMaker: virtual public Detector{
  public:
    static std::unique_ptr<DetectorSystem> make(const SimpleConfig&);
  };

} //namespace mu2e

#endif /* GeometryService_src_DetectorSystemMaker_hh */
