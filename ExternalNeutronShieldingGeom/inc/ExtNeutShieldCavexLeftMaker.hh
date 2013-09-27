#ifndef ExternalNeutronShieldingGeom_ExtNeutShieldCavexLeftMaker_hh
#define ExternalNeutronShieldingGeom_ExtNeutShieldCavexLeftMaker_hh
//
// Class to construct and return ExtNeutShieldCavexLeft
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class ExtNeutShieldCavexLeft;
  class SimpleConfig;

  class ExtNeutShieldCavexLeftMaker {
  public:

    static std::unique_ptr<ExtNeutShieldCavexLeft>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ExternalNeutronShieldingGeom_ExtNeutShieldCavexLeftMaker_hh */
