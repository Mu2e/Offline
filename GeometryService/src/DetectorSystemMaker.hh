#ifndef DetectorSystemMaker_HH
#define DetectorSystemMaker_HH
//
// Construct a DetectorSystem object.
//
// $Id: DetectorSystemMaker.hh,v 1.1 2010/12/03 00:52:40 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/12/03 00:52:40 $
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

#endif
