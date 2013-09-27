#ifndef ExternalNeutronShieldingGeom_ExtNeutShieldLAboveMaker_hh
#define ExternalNeutronShieldingGeom_ExtNeutShieldLAboveMaker_hh
//
// Class to construct and return ExtNeutShieldLAbove
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class ExtNeutShieldLAbove;
  class SimpleConfig;

  class ExtNeutShieldLAboveMaker {
  public:

    static std::unique_ptr<ExtNeutShieldLAbove>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ExternalNeutronShieldingGeom_ExtNeutShieldLAboveMaker_hh */
