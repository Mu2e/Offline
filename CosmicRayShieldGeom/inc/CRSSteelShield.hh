#ifndef CRSSteelShield_hh
#define CRSSteelShield_hh

//
// Representation of CRSSteelShield aka the flux return yoke
//
// $Id: CRSSteelShield.hh,v 1.1 2011/03/09 19:24:29 genser Exp $
// $Author: genser $ 
// $Date: 2011/03/09 19:24:29 $
//
// Original author KLG
//
// c++ includes
#include <vector>
#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class CRSSteelShield {

    friend class CosmicRayShieldMaker;

  public:

    CRSSteelShield() :
      _name(),
      _localOffset(),
      _rotation(),
      _globalOffset(),
      _halfLengths(),
      _holeRadius()
    {};

    // accept default d'tor etc...

    CRSSteelShield(std::string const name,
                               CLHEP::Hep3Vector   localOffset,
                               CLHEP::HepRotation* rotation,
                               CLHEP::Hep3Vector   globalOffset,
                               double const        halfLengths[3],
                               double              holeRadius = 0.):
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

    std::string         const & name() const { return _name;};
    CLHEP::Hep3Vector   const & getLocalOffset()  const { return _localOffset; };
    CLHEP::Hep3Vector   const & getGlobalOffset() const { return _globalOffset; };
    CLHEP::HepRotation* const  getRotation()     const { return _rotation; };

    // Return halfLengths, yes nestBox can take it directly as an argument
    double const* getHalfLengths()  const { return _halfLengths; };

    // we assume that the hole is in the center of the given shield...
    double const & getHoleRadius()                const { return _holeRadius; };

  private:

    std::string _name;

    // position of the center in the parent frame
    CLHEP::Hep3Vector _localOffset;

    // unit vector along the shield=bars direction
    CLHEP::Hep3Vector _uv; // not used for now

    // rotation in the parent frame; 0 for now
    CLHEP::HepRotation* _rotation;

    // position of the center in the global Mu2e frame
    CLHEP::Hep3Vector _globalOffset;

    double _halfLengths[3];

    double _holeRadius;

  };

}

#endif
