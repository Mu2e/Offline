#ifndef GeometryService_src_DetectorSystemMaker_hh
#define GeometryService_src_DetectorSystemMaker_hh
//
// Construct a DetectorSystem object.
//
// $Id: DetectorSystemMaker.hh,v 1.3 2012/02/24 16:36:36 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/24 16:36:36 $
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

  class DetectorSystemMaker: public Detector{

  public:

    DetectorSystemMaker( SimpleConfig const& );

    // Accept the compiler generator d'tor, copy c'tor and copy assigment.

    // Give control of the detector system to the caller.
    std::auto_ptr<DetectorSystem> getDetectorSystemPtr() { return _detectorSystem; }

  private:

    std::auto_ptr<DetectorSystem> _detectorSystem;

  };

} //namespace mu2e

#endif /* GeometryService_src_DetectorSystemMaker_hh */
