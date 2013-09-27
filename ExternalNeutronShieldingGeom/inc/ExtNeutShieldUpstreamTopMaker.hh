#ifndef ExternalNeutronShieldingGeom_ExtNeutShieldUpstreamTopMaker_hh
#define ExternalNeutronShieldingGeom_ExtNeutShieldUpstreamTopMaker_hh
//
// Class to construct and return ExtNeutShieldUpstreamTop
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class ExtNeutShieldUpstreamTop;
  class SimpleConfig;

  class ExtNeutShieldUpstreamTopMaker {
  public:

    static std::unique_ptr<ExtNeutShieldUpstreamTop>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ExternalNeutronShieldingGeom_ExtNeutShieldUpstreamTopMaker_hh */
