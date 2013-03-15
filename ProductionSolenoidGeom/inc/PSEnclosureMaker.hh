#ifndef ProductionSolenoidGeom_PSEnclosureMaker_hh
#define ProductionSolenoidGeom_PSEnclosureMaker_hh
//
// Class to construct and return PSEnclosure
//
// $Id: PSEnclosureMaker.hh,v 1.5 2013/03/15 15:52:05 kutschke Exp $
// $Author: kutschke $
// $Date: 2013/03/15 15:52:05 $
//
// Original author Andrei Gaponenko
//

#include <string>
#include <memory>

namespace CLHEP { class Hep3Vector; }

namespace mu2e {

  class PSEnclosure;
  class SimpleConfig;

  class PSEnclosureMaker {
  public:

    static std::unique_ptr<PSEnclosure>  make(const SimpleConfig& config,

                                            // The center of the downstream surface of the PS
                                            const CLHEP::Hep3Vector& psEndRefPoint);
  };

}  //namespace mu2e

#endif /* ProductionSolenoidGeom_PSEnclosureMaker_hh */
