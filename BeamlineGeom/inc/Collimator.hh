#ifndef BeamlineGeom_Collimator_hh
#define BeamlineGeom_Collimator_hh

//
// Class to represent the collimators. Position is relative to
// the corresponding TS straight section
//
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class Collimator {

  friend class BeamlineMaker;

  public:

    Collimator() : _halfZ(0.0), _origin() {}

    // use compiler-generated copy c'tor, copy assignment, and d'tor

    Collimator(double halfZ, CLHEP::Hep3Vector const& origin) :  _halfZ(halfZ), _origin(origin) {}

    void set(double halfZ, CLHEP::Hep3Vector const& origin) {
      _halfZ=halfZ;
      _origin=origin;
    }

    double halfLength() const { return _halfZ; }
    CLHEP::Hep3Vector const& getLocal() const { return _origin; }

  private:

    double _halfZ;
    CLHEP::Hep3Vector _origin;

    void adjustOffset( const CLHEP::Hep3Vector& offset ) {
      _origin=_origin+offset;
    }

};

}
#endif /* BeamlineGeom_Collimator_hh */
