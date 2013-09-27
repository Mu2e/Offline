#ifndef ExternalNeutronShieldingGeom_ExtNeutShieldLCeilingMaker_hh
#define ExternalNeutronShieldingGeom_ExtNeutShieldLCeilingMaker_hh
//
// Class to construct and return ExtNeutShieldLCeiling
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class ExtNeutShieldLCeiling;
  class SimpleConfig;

  class ExtNeutShieldLCeilingMaker {
  public:

    static std::unique_ptr<ExtNeutShieldLCeiling>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ExternalNeutronShieldingGeom_ExtNeutShieldLCeilingMaker_hh */
