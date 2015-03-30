#ifndef ExternalShieldingGeom_SaddleMaker_hh
#define ExternalShieldingGeom_SaddleMaker_hh
//
// Class to construct and return Saddle
//
//
// Original author David Norvil Brown
//

#include <string>
#include <memory>


namespace mu2e {

  class Saddle;
  class SimpleConfig;

  class SaddleMaker {
  public:

    static std::unique_ptr<Saddle>  make(const SimpleConfig& 
						      config );

  };

}  //namespace mu2e

#endif /* ExternalShieldingGeom_SaddleMaker_hh */
