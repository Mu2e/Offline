#ifndef ProductionSolenoidGeom_PSExternalShieldingMaker_hh
#define ProductionSolenoidGeom_PSExternalShieldingMaker_hh
//
// Class to construct and return PSExternalShielding
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class PSExternalShielding;
  class SimpleConfig;

  class PSExternalShieldingMaker {
  public:

    static std::unique_ptr<PSExternalShielding>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ProductionSolenoidGeom_PSExternalShieldingMaker_hh */
