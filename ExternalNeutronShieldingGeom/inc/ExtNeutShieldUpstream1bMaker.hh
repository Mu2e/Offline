#ifndef ExternalNeutronShieldingGeom_ExtNeutShieldUpstream1bMaker_hh
#define ExternalNeutronShieldingGeom_ExtNeutShieldUpstream1bMaker_hh
//
// Class to construct and return ExtNeutShieldUpstream1b
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class ExtNeutShieldUpstream1b;
  class SimpleConfig;

  class ExtNeutShieldUpstream1bMaker {
  public:

    static std::unique_ptr<ExtNeutShieldUpstream1b>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ExternalNeutronShieldingGeom_ExtNeutShieldUpstream1bMaker_hh */
