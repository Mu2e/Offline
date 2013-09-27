#ifndef ExternalNeutronShieldingGeom_ExtNeutShieldCendBoxesMaker_hh
#define ExternalNeutronShieldingGeom_ExtNeutShieldCendBoxesMaker_hh
//
// Class to construct and return ExtNeutShieldCendBoxes
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class ExtNeutShieldCendBoxes;
  class SimpleConfig;

  class ExtNeutShieldCendBoxesMaker {
  public:

    static std::unique_ptr<ExtNeutShieldCendBoxes>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ExternalNeutronShieldingGeom_ExtNeutShieldCendBoxesMaker_hh */
