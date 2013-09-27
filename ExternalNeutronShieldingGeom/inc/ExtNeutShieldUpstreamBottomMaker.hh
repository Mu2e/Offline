#ifndef ExternalNeutronShieldingGeom_ExtNeutShieldUpstreamBottomMaker_hh
#define ExternalNeutronShieldingGeom_ExtNeutShieldUpstreamBottomMaker_hh
//
// Class to construct and return ExtNeutShieldUpstreamBottom
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class ExtNeutShieldUpstreamBottom;
  class SimpleConfig;

  class ExtNeutShieldUpstreamBottomMaker {
  public:

    static std::unique_ptr<ExtNeutShieldUpstreamBottom>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ExternalNeutronShieldingGeom_ExtNeutShieldUpstreamBottomMaker_hh */
