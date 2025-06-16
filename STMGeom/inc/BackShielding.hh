#ifndef STMGeom_BackShielding_hh
#define STMGeom_BackShielding_hh

// Germanium Detector Object
//
// Author: Haichuan Cao
// Sept 2023

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class BackShielding {
  public:
    BackShielding(bool build, double thick, double length, double height,
    double Back_dX, double Back_dY, double STMShieldingPipeGap,
    CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(), CLHEP::HepRotation const & rotation = CLHEP::HepRotation()
    ):
      _build(build), _BPThick(thick), _BPLength(length), _BPHeight(height),
      _Back_dX(Back_dX), _Back_dY(Back_dY), _STMShieldingPipeGap(STMShieldingPipeGap),
      _originInMu2e(originInMu2e), _rotation(rotation)
    {}

   bool     build()                    const {return _build;}
   double   BPThick()                  const {return _BPThick;}
   double   BPLength()                 const {return _BPLength;}
   double   BPHeight()                 const {return _BPHeight;}
   double   Back_dX()                  const {return _Back_dX;}
   double   Back_dY()                  const {return _Back_dY;}
   double   STMShieldingPipeGap()                  const {return _STMShieldingPipeGap;}

    BackShielding() {}
  private:

    bool               _build;
    double             _BPThick;
    double             _BPLength;
    double             _BPHeight;
    double             _Back_dX;
    double             _Back_dY;
    double             _STMShieldingPipeGap;

    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation;
  };

}

#endif/*STMGeom_BackShielding_hh*/
