#ifndef ExternalNeutronShieldingGeom_ExtNeutShieldUpstream2Maker_hh
#define ExternalNeutronShieldingGeom_ExtNeutShieldUpstream2Maker_hh
//
// Class to construct and return ExtNeutShieldUpstream2
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class ExtNeutShieldUpstream2;
  class SimpleConfig;

  class ExtNeutShieldUpstream2Maker {
  public:

    static std::unique_ptr<ExtNeutShieldUpstream2>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ExternalNeutronShieldingGeom_ExtNeutShieldUpstream2Maker_hh */
