#ifndef ProductionSolenoidGeom_PSShieldMaker_hh
#define ProductionSolenoidGeom_PSShieldMaker_hh
//
// Class to construct and return PSShield
//
// $Id: PSShieldMaker.hh,v 1.2 2012/04/27 05:37:32 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/27 05:37:32 $
//
// Original author Andrei Gaponenko
//

#include <memory>
#include "ProductionSolenoidGeom/inc/PSShield.hh"

namespace CLHEP { class Hep3Vector; }

namespace mu2e {

  class SimpleConfig;

  class PSShieldMaker {
    static PSShield::Groove readGroove(int i, const SimpleConfig& config);

  public:

    static std::auto_ptr<PSShield> make(const SimpleConfig& config,
                                        // The center of the downstream surface of the PS
                                        const CLHEP::Hep3Vector& psEndRefPoint);
  };

}  //namespace mu2e

#endif /* ProductionSolenoidGeom_PSShieldMaker_hh */
