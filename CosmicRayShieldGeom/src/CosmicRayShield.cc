//
// Class to represent CosmicRayShield
//
// $Id: CosmicRayShield.cc,v 1.1 2011/01/25 16:43:52 genser Exp $
// $Author: genser $ 
// $Date: 2011/01/25 16:43:52 $
//
// Original author KLG
//

#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"

namespace mu2e {

  // Get SteelShield
  CosmicRayShieldSteelShield const & CosmicRayShield::getCosmicRayShieldSteelShield(std::string name) const { 
    return _steelShield.find(name)->second;
  };
  
  // Get Veto
  // CosmicRayShieldVeto  const& getCosmicRayShieldVeto()  const { return _veto; };

  void CosmicRayShield::addSteelShield(std::string name, 
                                       CLHEP::Hep3Vector   localOffset, 
                                       CLHEP::HepRotation* localRot,
                                       CLHEP::Hep3Vector   globalOffset,
                                       double const        halfLengths[3],
                                       double              holeRadius) 
  {
    
    _steelShield[name] = CosmicRayShieldSteelShield(name,
                                                    localOffset,
                                                    localRot,
                                                    globalOffset,
                                                    halfLengths,
                                                    holeRadius);
  };

}

