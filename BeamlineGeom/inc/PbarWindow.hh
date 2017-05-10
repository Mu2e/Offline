#ifndef BeamlineGeom_PbarWindow_hh
#define BeamlineGeom_PbarWindow_hh

//
// Class to represent the collimators. Position is relative to
// the corresponding TS straight section
//
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class PbarWindow {

  friend class BeamlineMaker;

  public:

    PbarWindow() : _rOut(0.), _halfZ(0.) {}

    // use compiler-generated copy c'tor, copy assignment, and d'tor

    std::string shape() const { return _shape; }
    double halfLength() const { return _halfZ; }
    CLHEP::Hep3Vector const& getLocal() const { return _origin; }
    std::string material() const { return _material; }

    double rOut()   const  { return _rOut; }
    double getY0()  const  { return _y0 ; };
    double getY1()  const  { return _y1 ; };
    double getDZ0() const  { return _dz0; };
    double getDZ1() const  { return _dz1; };
    double getWedgeZOffset() const { return _wedgeZOffset; }

    void set(double halfZ, CLHEP::Hep3Vector origin) {
      _halfZ  = halfZ;
      _origin = origin; 
    }

  private:

    std::string _shape;
    std::string _material;

    double _rOut;
    double _halfZ;
    CLHEP::Hep3Vector _origin;

    double _y0 ;
    double _y1 ;
    double _dz0;
    double _dz1;
    double _wedgeZOffset;
};

}
#endif /* BeamlineGeom_PbarWindow_hh */
