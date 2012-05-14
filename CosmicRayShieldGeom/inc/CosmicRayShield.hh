#ifndef CosmicRayShieldGeom_CosmicRayShield_hh
#define CosmicRayShieldGeom_CosmicRayShield_hh

//
// Representation of CosmicRayShield
//
// $Id: CosmicRayShield.hh,v 1.14 2012/05/14 21:22:24 genser Exp $
// $Author: genser $
// $Date: 2012/05/14 21:22:24 $
//
// Original author KLG
//

// c++ includes
#include <map>
#include <string>

// clhep includes
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

// Includes from Mu2e
#include "CosmicRayShieldGeom/inc/CRSScintillatorShield.hh"
#include "CosmicRayShieldGeom/inc/CRSSteelShield.hh"
#include "Mu2eInterfaces/inc/Detector.hh"


namespace mu2e {

  // Forward reference.
  class SimpleConfig;
  class CosmicRayShieldMaker;

  class CosmicRayShield : virtual public Detector {

    friend class CosmicRayShieldMaker;

  public:

    CosmicRayShield() {
      _hasActiveShield = false;
      _hasPassiveShield = false;
    }

    ~CosmicRayShield(){}

    // Get SteelShield
    CRSSteelShield         const & getCRSSteelShield(std::string name) const;

    std::map<std::string,CRSSteelShield> const & getCRSSteelShields() const {
      return _steelShields;
    }

    // Get ScintillatorShield
    CRSScintillatorShield  const & getCRSScintillatorShield(std::string name)  const;

    std::map<std::string,CRSScintillatorShield> const & getCRSScintillatorShields() const {
      return _scintillatorShields;
    }

    CRSScintillatorBarDetail const & getCRSScintillatorBarDetail() const {
      return _barDetails;
    }

    std::vector<CRSScintillatorBar> const & getAllCRSScintillatorBars() const {
      return _allCRSScintillatorBars;
    }

    const CRSScintillatorBar& getBar ( CRSScintillatorBarIndex index ) const {
      return _allCRSScintillatorBars.at(index.asInt());
    }

    bool const hasActiveShield() const {
      return _hasActiveShield;
    }

    bool const hasPassiveShield() const {
      return _hasPassiveShield;
    }

  private:

    std::map<std::string,CRSSteelShield>         _steelShields;

    std::map<std::string,CRSScintillatorShield>  _scintillatorShields;

    // Detailed info about scintillators etc...
    CRSScintillatorBarDetail _barDetails;

    // global holder of all scintillator bars
    std::vector<CRSScintillatorBar>  _allCRSScintillatorBars;

    bool _hasActiveShield;
    bool _hasPassiveShield;

    // for a future reference:

    // given that the bars, layers and modules have the same rotation
    // angles as their shield, one could e.g. store the rotations as
    // HepRotation in the shields or in CosmicRayShield and point to
    // them from the lower level objects; similarly for their half
    // lengths, one could store layer and module half lengths as "details"
    // in CosmicRayShield and point to them from those objects

  };

}
#endif /* CosmicRayShieldGeom_CosmicRayShield_hh */
