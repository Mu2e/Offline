#ifndef GeometryService_src_DetectorSystemMaker_hh
#define GeometryService_src_DetectorSystemMaker_hh
//
// Construct a DetectorSystem object.
//
// $Id: DetectorSystemMaker.hh,v 1.6 2013/03/15 15:52:04 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:04 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <memory>

// Mu2e includes.
#include "Mu2eInterfaces/inc/Detector.hh"
#include "GeometryService/inc/DetectorSystem.hh"

namespace mu2e {

  // Forward references
  class SimpleConfig;

  class DetectorSystemMaker: virtual public Detector{
  public:
    static std::unique_ptr<DetectorSystem> make(const SimpleConfig&);
  };

} //namespace mu2e

#endif /* GeometryService_src_DetectorSystemMaker_hh */
