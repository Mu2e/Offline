#ifndef ExternalNeutronShieldingGeom_ExtNeutShieldCryoBoxesMaker_hh
#define ExternalNeutronShieldingGeom_ExtNeutShieldCryoBoxesMaker_hh
//
// Class to construct and return ExtNeutShieldCryoBoxes
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class ExtNeutShieldCryoBoxes;
  class SimpleConfig;

  class ExtNeutShieldCryoBoxesMaker {
  public:

    static std::unique_ptr<ExtNeutShieldCryoBoxes>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ExternalNeutronShieldingGeom_ExtNeutShieldCryoBoxesMaker_hh */
