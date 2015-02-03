#ifndef ExternalShieldingGeom_ExtShieldDownstreamMaker_hh
#define ExternalShieldingGeom_ExtShieldDownstreamMaker_hh
//
// Class to construct and return ExtShieldDownstream
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class ExtShieldDownstream;
  class SimpleConfig;

  class ExtShieldDownstreamMaker {
  public:

    static std::unique_ptr<ExtShieldDownstream>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ExternalShieldingGeom_ExtShieldDownstreamMaker_hh */
