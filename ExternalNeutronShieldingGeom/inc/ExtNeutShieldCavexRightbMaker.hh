#ifndef ExternalNeutronShieldingGeom_ExtNeutShieldCavexRightbMaker_hh
#define ExternalNeutronShieldingGeom_ExtNeutShieldCavexRightbMaker_hh
//
// Class to construct and return ExtNeutShieldCavexRightb
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class ExtNeutShieldCavexRightb;
  class SimpleConfig;

  class ExtNeutShieldCavexRightbMaker {
  public:

    static std::unique_ptr<ExtNeutShieldCavexRightb>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ExternalNeutronShieldingGeom_ExtNeutShieldCavexRightbMaker_hh */
