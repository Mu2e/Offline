#ifndef Collimator_HH
#define Collimator_HH

//
// Class to represent the collimators. Position is relative to
// the corresponding TS straight section
//
#include <string>
#include <memory>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class Collimator {

  friend class BeamlineMaker;

  public:

    Collimator() { ; };

    explicit Collimator(double halfZ, CLHEP::Hep3Vector origin) :  _halfZ(halfZ), _origin(origin) {;};

    ~Collimator(){ };

    void set(double halfZ, CLHEP::Hep3Vector origin) {
      _halfZ=halfZ;
      _origin=origin;
    };

    double getHalfLength() const { return _halfZ; }
    CLHEP::Hep3Vector const& getLocal() const { return _origin; }

  private:

    double _halfZ;
    CLHEP::Hep3Vector _origin;

};

}
#endif
