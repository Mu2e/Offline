#ifndef ProductionSolenoidGeom_PSEnclosureMaker_hh
#define ProductionSolenoidGeom_PSEnclosureMaker_hh
//
// Class to construct and return PSEnclosure
//
// $Id: PSEnclosureMaker.hh,v 1.3 2012/03/30 19:18:03 gandr Exp $
// $Author: gandr $
// $Date: 2012/03/30 19:18:03 $
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

    static std::auto_ptr<PSEnclosure>  make(const SimpleConfig& config,

                                            // The center of the downstream surface of the PS
                                            const CLHEP::Hep3Vector& psEndRefPoint,

                                            const std::string& psInsideMaterialName);
  };

}  //namespace mu2e

#endif /* ProductionSolenoidGeom_PSEnclosureMaker_hh */
