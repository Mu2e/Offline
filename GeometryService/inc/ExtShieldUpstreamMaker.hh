#ifndef ExternalShieldingGeom_ExtShieldUpstreamMaker_hh
#define ExternalShieldingGeom_ExtShieldUpstreamMaker_hh
//
// Class to construct and return ExtShieldUpstream
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class ExtShieldUpstream;
  class SimpleConfig;

  class ExtShieldUpstreamMaker {
  public:

    static std::unique_ptr<ExtShieldUpstream>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ExternalShieldingGeom_ExtShieldUpstreamMaker_hh */
