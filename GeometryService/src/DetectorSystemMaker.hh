#ifndef GeometryService_src_DetectorSystemMaker_hh
#define GeometryService_src_DetectorSystemMaker_hh
//
// Construct a DetectorSystem object.
//
// $Id: DetectorSystemMaker.hh,v 1.5 2012/03/30 19:18:03 gandr Exp $
// $Author: gandr $
// $Date: 2012/03/30 19:18:03 $
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
    static std::auto_ptr<DetectorSystem> make(const SimpleConfig&);
  };

} //namespace mu2e

#endif /* GeometryService_src_DetectorSystemMaker_hh */
