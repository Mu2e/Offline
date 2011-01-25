#ifndef CosmicRayShield_hh
#define CosmicRayShield_hh

//
// Class to represent CosmicRayShield
//
// $Id: CosmicRayShield.hh,v 1.1 2011/01/25 16:43:52 genser Exp $
// $Author: genser $ 
// $Date: 2011/01/25 16:43:52 $
//
// Original author KLG
//

// c++ includes
#include <string>
#include <map>

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

// Includes from Mu2e
#include "GeometryService/inc/Detector.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShieldSteelShield.hh"


namespace mu2e {

  // Forward reference.
  class SimpleConfig;
  class CosmicRayShieldMaker;

  class CosmicRayShield : public Detector {

    friend class CosmicRayShieldMaker;
    //  friend class CosmicRayShieldVetoMaker;

  public:

    CosmicRayShield():
      _name("CosmicRayShield")
    {};

    ~CosmicRayShield(){};

    virtual std::string name() const { return _name;};
    
    // Get SteelShield
    CosmicRayShieldSteelShield const & getCosmicRayShieldSteelShield(std::string name) const;

    // Get Veto
    //    CosmicRayShieldVeto  const& getCosmicRayShieldVeto()  const { return _veto; };

    void addSteelShield(std::string name, 
                        CLHEP::Hep3Vector   localOffset, 
                        CLHEP::HepRotation* localRot,
                        CLHEP::Hep3Vector   globalOffset,
                        double const        halfLengths[3],
                        double              holeRadius = 0.);

  private:

    std::string _name;

    std::map<std::string,CosmicRayShieldSteelShield> _steelShield;

    //    CosmicRayShieldVeto  _veto;

  };

}
#endif
