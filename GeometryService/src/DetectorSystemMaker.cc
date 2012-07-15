//
// Construct a DetectorSystem object.
//
// $Id: DetectorSystemMaker.cc,v 1.3 2012/07/15 22:06:17 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/15 22:06:17 $
//
// Original author Rob Kutschke
//

// Mu2e includes.
#include "GeometryService/src/DetectorSystemMaker.hh"
#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e{

  std::auto_ptr<DetectorSystem> DetectorSystemMaker::make(const SimpleConfig& config) {

    // The detector system origin, as measured in the Mu2e system, is on the
    // axis of the DS and is at the specified z.
    CLHEP::Hep3Vector origin( - config.getDouble("mu2e.solenoidOffset"),
                                0.,
                                config.getDouble("mu2e.detectorSystemZ0")
                              );

    return std::auto_ptr<DetectorSystem>(new DetectorSystem(origin));
  }

}
