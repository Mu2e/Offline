#ifndef CosmicRayShieldSteelShield_hh
#define CosmicRayShieldSteelShield_hh

//
// Class to represent CosmicRayShieldSteelShield
//
// $Id: CosmicRayShieldSteelShield.hh,v 1.1 2011/01/25 16:43:52 genser Exp $
// $Author: genser $ 
// $Date: 2011/01/25 16:43:52 $
//
// Original author KLG
//
// c++ includes
#include <vector>
#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  // Forward reference.
  class SimpleConfig;

  class CosmicRayShieldSteelShield {

    friend class CosmicRayShieldSteelMaker;

  public:

    CosmicRayShieldSteelShield() :
      _name(),
      _localOffset(),
      _rotation(),
      _globalOffset(),
      _halfLengths(),
      _holeRadius()
    {};

    // accept default d'tor etc...

    // CLHEP::HepRotation should be "owned" somwhere else, they are "0' for now anyway...

    CosmicRayShieldSteelShield(std::string const name,
                               CLHEP::Hep3Vector localOffset,
                               CLHEP::HepRotation* rotation,
                               CLHEP::Hep3Vector  globalOffset,
                               double const halfLengths[3],
                               double holeRadius):
      _name(name),
      _localOffset(localOffset),
      _rotation(rotation),
      _globalOffset(globalOffset),
      _holeRadius(holeRadius)
    {      
      _halfLengths[0] = halfLengths[0];
      _halfLengths[1] = halfLengths[1];
      _halfLengths[2] = halfLengths[2];
    };

    std::string         const& name() const { return _name;};
    CLHEP::Hep3Vector   const& getLocalOffset()  const { return _localOffset; };
    CLHEP::Hep3Vector   const& getGlobalOffset() const { return _globalOffset; };
    CLHEP::HepRotation* const  getRotation()     const { return _rotation; };

    // Return halfLengths, yes nestBox can take it directly as an argument
    double const* getHalfLengths()  const { return _halfLengths; };

    // we assume that the hole is in the center of the given shield...
    double const& getHoleRadius()                const { return _holeRadius; };

  private:

    // one could also create an object which would contain one of each
    // below..., e.g. contained in the object above it

    std::string _name;

    // position in the parent frame
    CLHEP::Hep3Vector _localOffset;

    // rotation in the parent frame
    CLHEP::HepRotation* _rotation;

    // position in the global Mu2e frame
    CLHEP::Hep3Vector _globalOffset;

    double _halfLengths[3];

    double _holeRadius;

  };

}

#endif
