#ifndef GeometryService_src_DetectorSystemMaker_hh
#define GeometryService_src_DetectorSystemMaker_hh
//
// Construct a DetectorSystem object.
//
// $Id: DetectorSystemMaker.hh,v 1.2 2011/05/17 15:41:35 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:35 $
//
// Original author Rob Kutschke
//

// C++ includes
#include <memory>

// Mu2e includes.
#include "GeometryService/inc/Detector.hh"
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
