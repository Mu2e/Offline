#ifndef GeomPrimitives_Torus_hh
#define GeomPrimitives_Torus_hh
//
// "Typical" Torus object
//
//
// Original author: Kyle Knoepfel
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e {

  class Torus {

  public:

    Torus(double torusRadius, double innerRadius, double outerRadius, 
          double phi0, double deltaPhi,
          CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(),
          CLHEP::HepRotation const & rotation = CLHEP::HepRotation(), 
          std::string const & materialName = "" ) :
      _torusRadius( torusRadius ),
      _innerRadius( innerRadius ),
      _outerRadius( outerRadius ),
      _phi0       ( phi0        ),
      _deltaPhi   ( deltaPhi    ),
      _originInMu2e( originInMu2e ),
      _rotation    ( rotation     ),
      _materialName( materialName )
    {}

    double torusRadius() const { return _torusRadius; }
    double innerRadius() const { return _innerRadius; }
    double outerRadius() const { return _outerRadius; }
    double phi0()        const { return _phi0       ; }
    double deltaPhi()    const { return _deltaPhi   ; }

    std::string const & materialName() const { return _materialName; }

    CLHEP::Hep3Vector const  & originInMu2e() const { return _originInMu2e; }
    CLHEP::HepRotation const & rotation()     const { return _rotation; }

    // Genreflex can't do persistency of vector<Torus> without a default constructor
    Torus() {}

  private:

    double _torusRadius; 
    double _innerRadius; 
    double _outerRadius; 
    double _phi0       ; 
    double _deltaPhi   ; 

    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation; // wrt to parent volume
    std::string        _materialName;

  };

}

#endif/*GeomPrimitives_Torus_hh*/
