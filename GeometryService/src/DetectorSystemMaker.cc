//
// Construct a DetectorSystem object.
//
// $Id: DetectorSystemMaker.cc,v 1.1 2010/12/03 00:52:40 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/12/03 00:52:40 $
//
// Original author Rob Kutschke
//

// Mu2e includes.
#include "GeometryService/src/DetectorSystemMaker.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

namespace mu2e{

  DetectorSystemMaker::DetectorSystemMaker( SimpleConfig const& config ){


    // The detector system origin, as measured in the Mu2e system, is on the
    // axis of the DS and is at the specified z.
    CLHEP::Hep3Vector origin( - config.getDouble("mu2e.solenoidOffset"),
                                0.,
                                config.getDouble("mu2e.detectorSystemZ0")
                              );

    _detectorSystem = std::auto_ptr<DetectorSystem>( new DetectorSystem(origin) );

  }

}
