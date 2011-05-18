//
// Class to represent CosmicRayShield
//
// $Id: CosmicRayShield.cc,v 1.3 2011/05/18 02:27:15 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:15 $
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


