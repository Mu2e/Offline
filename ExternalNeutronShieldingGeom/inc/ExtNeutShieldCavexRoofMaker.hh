#ifndef ExternalNeutronShieldingGeom_ExtNeutShieldCavexRoofMaker_hh
#define ExternalNeutronShieldingGeom_ExtNeutShieldCavexRoofMaker_hh
//
// Class to construct and return ExtNeutShieldCavexRoof
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class ExtNeutShieldCavexRoof;
  class SimpleConfig;

  class ExtNeutShieldCavexRoofMaker {
  public:

    static std::unique_ptr<ExtNeutShieldCavexRoof>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ExternalNeutronShieldingGeom_ExtNeutShieldCavexRoofMaker_hh */
