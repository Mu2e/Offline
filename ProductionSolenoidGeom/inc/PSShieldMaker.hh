#ifndef ProductionSolenoidGeom_PSShieldMaker_hh
#define ProductionSolenoidGeom_PSShieldMaker_hh
//
// Class to construct and return PSShield
//
// $Id: PSShieldMaker.hh,v 1.1 2012/03/30 14:07:14 gandr Exp $
// $Author: gandr $
// $Date: 2012/03/30 14:07:14 $
//
// Original author Andrei Gaponenko
//

#include <memory>

namespace CLHEP { class Hep3Vector; }

namespace mu2e {

  class PSShield;
  class SimpleConfig;

  class PSShieldMaker {
  public:

    static std::auto_ptr<PSShield> make(const SimpleConfig& config,
                                        // The center of the downstream surface of the PS
                                        const CLHEP::Hep3Vector& psEndRefPoint);
  };

}  //namespace mu2e

#endif /* ProductionSolenoidGeom_PSShieldMaker_hh */
