//
// Class to represent CosmicRayShield
//
// $Id: CosmicRayShield.cc,v 1.4 2013/09/13 06:42:44 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2013/09/13 06:42:44 $
//
// Original author KLG
//

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"

namespace mu2e 
{
  // Get ScintillatorShield
  CRSScintillatorShield const &CosmicRayShield::getCRSScintillatorShield(std::string name) const 
  {
    return _scintillatorShields.find(name)->second;
  }
}


