#ifndef ExternalNeutronShieldingGeom_ExtNeutShieldUpstream1aMaker_hh
#define ExternalNeutronShieldingGeom_ExtNeutShieldUpstream1aMaker_hh
//
// Class to construct and return ExtNeutShieldUpstream1a
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class ExtNeutShieldUpstream1a;
  class SimpleConfig;

  class ExtNeutShieldUpstream1aMaker {
  public:

    static std::unique_ptr<ExtNeutShieldUpstream1a>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ExternalNeutronShieldingGeom_ExtNeutShieldUpstream1aMaker_hh */
