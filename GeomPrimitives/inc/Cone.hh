#ifndef GeomPrimitives_Cone_hh
#define GeomPrimitives_Cone_hh
//
// "Typical" Cone object
//
// Original author: KLG
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
//#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e {

  class Cone {

  public:

    Cone(double innerRadius1, double outerRadius1,
         double innerRadius2, double outerRadius2,
         double halfLength,
         double phi0, double deltaPhi,
         CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(),
         CLHEP::HepRotation const & rotation = CLHEP::HepRotation(), 
         std::string const & materialName = "" ) :
      _innerRadius1( innerRadius1 ),
      _outerRadius1( outerRadius1 ),
      _innerRadius2( innerRadius2 ),
      _outerRadius2( outerRadius2 ),
      _halfLength  ( halfLength   ),
      _phi0        ( phi0         ),
      _deltaPhi    ( deltaPhi     ),
      _originInMu2e( originInMu2e ),
      _rotation    ( rotation     ),
      _materialName( materialName )
    {}

    double innerRadius1() const { return _innerRadius1; }
    double outerRadius1() const { return _outerRadius1; }
    double innerRadius2() const { return _innerRadius2; }
    double outerRadius2() const { return _outerRadius2; }
    double halfLength()   const { return _halfLength; }
    double zHalfLength()  const { return _halfLength; }
    double zBegin()       const { return _originInMu2e.z() - zHalfLength(); }
    double zEnd()         const { return _originInMu2e.z() + zHalfLength(); }
    double phi0()         const { return _phi0; }
    double deltaPhi()     const { return _deltaPhi; }

    std::string const & materialName() const { return _materialName; }

    CLHEP::Hep3Vector const  & originInMu2e() const { return _originInMu2e; }
    CLHEP::HepRotation const & rotation()     const { return _rotation; }

    // Genreflex can't do persistency of vector<Cone> without a default constructor
    Cone() {}

  private:

    double _innerRadius1; // at -pDz
    double _outerRadius1;
    double _innerRadius2; // at +pDz
    double _outerRadius2;
    double _halfLength;
    double _phi0; 
    double _deltaPhi; 

    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation; // wrt to parent volume
    std::string        _materialName;

  };

}

#endif/*GeomPrimitives_Cone_hh*/
