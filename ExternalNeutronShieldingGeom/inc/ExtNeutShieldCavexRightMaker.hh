#ifndef ExternalNeutronShieldingGeom_ExtNeutShieldCavexRightMaker_hh
#define ExternalNeutronShieldingGeom_ExtNeutShieldCavexRightMaker_hh
//
// Class to construct and return ExtNeutShieldCavexRight
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class ExtNeutShieldCavexRight;
  class SimpleConfig;

  class ExtNeutShieldCavexRightMaker {
  public:

    static std::unique_ptr<ExtNeutShieldCavexRight>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ExternalNeutronShieldingGeom_ExtNeutShieldCavexRightMaker_hh */
