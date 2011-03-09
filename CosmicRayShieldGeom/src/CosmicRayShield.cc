//
// Class to represent CosmicRayShield
//
// $Id: CosmicRayShield.cc,v 1.2 2011/03/09 19:46:39 genser Exp $
// $Author: genser $ 
// $Date: 2011/03/09 19:46:39 $
//
// Original author KLG
//

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"

namespace mu2e {

  // Get SteelShield
  CRSSteelShield const & 
  CosmicRayShield::getCRSSteelShield(std::string name) const { 
    return _steelShields.find(name)->second;
  }
  
  // Get ScintillatorShield
  CRSScintillatorShield const & 
  CosmicRayShield::getCRSScintillatorShield(std::string name) const { 
    return _scintillatorShields.find(name)->second;
  }

}


