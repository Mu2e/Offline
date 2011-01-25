#ifndef CosmicRayShieldMaker_hh
#define CosmicRayShieldMaker_hh
//
// Construct and return CosmicRayShield.
//
// $Id: CosmicRayShieldMaker.hh,v 1.1 2011/01/25 16:43:52 genser Exp $
// $Author: genser $ 
// $Date: 2011/01/25 16:43:52 $
//
// Original author KLG
//

#include <vector>
#include <string>
#include <memory>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

class CosmicRayShield;
class SimpleConfig;

class CosmicRayShieldMaker {

public:

  CosmicRayShieldMaker( SimpleConfig const& config );  

  ~CosmicRayShieldMaker ();

  // This is depracted and will go away soon.  
  // Still needed for root graphics version.
  const CosmicRayShield& getCosmicRayShield() const { return *_crs;}

  // This is the accessor that will remain.
  std::auto_ptr<CosmicRayShield> getCosmicRayShieldPtr() { return _crs; }

private:

  std::auto_ptr<CosmicRayShield> _crs;

};

}  //namespace mu2e

#endif 
