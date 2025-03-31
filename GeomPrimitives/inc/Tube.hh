#ifndef GeomPrimitives_Tube_hh
#define GeomPrimitives_Tube_hh
//
// "Typical" Tube object
//
//
// Original author KLG
//

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "Offline/GeomPrimitives/inc/TubsParams.hh"


namespace mu2e {

  class TubsParams;
  class Tube {

  public:

    Tube(double rIn,double rOut,double halfLength, double phi0, double phiMax,
         std::string const & materialName, CLHEP::Hep3Vector const & originInMu2e,
         CLHEP::HepRotation const & rotation = CLHEP::HepRotation());

    Tube(double rIn,double rOut,double halfLength,
         CLHEP::Hep3Vector const & originInMu2e,
         CLHEP::HepRotation const & rotation = CLHEP::HepRotation(),
         double phi0 = 0., double phiMax = CLHEP::twopi,
         std::string const & materialName = "");

    Tube(std::string const & materialName, CLHEP::Hep3Vector const & originInMu2e,
         double rIn,double rOut,double halfLength,
         double phi0 = 0., double phiMax = CLHEP::twopi,
         CLHEP::HepRotation const & rotation = CLHEP::HepRotation() );

    double innerRadius() const { return _params.innerRadius(); }
    double outerRadius() const { return _params.outerRadius(); }
    double halfLength()  const { return _params.zHalfLength(); }
    double zHalfLength() const { return _params.zHalfLength(); }
    double zBegin()      const { return _originInMu2e.z() - zHalfLength() ; }
    double zEnd()        const { return _originInMu2e.z() + zHalfLength() ; }
    double phi0()        const { return _params.phi0(); }
    double phiMax()      const { return _params.phiMax(); }

    std::string const & materialName() const { return _materialName; }

    CLHEP::Hep3Vector const  & originInMu2e() const { return _originInMu2e; }
    CLHEP::HepRotation const & rotation()     const { return _rotation; }

    TubsParams const & getTubsParams() const { return _params; }

    std::string getName() const { return _name; }
    void        setName( std::string const& name ) { _name = name; }

    // Genreflex can't do persistency of vector<Tube> without a default constructor
    Tube() : _params(0,0,0) {}

  private:

    TubsParams         _params;
    CLHEP::Hep3Vector  _originInMu2e;
    std::string        _materialName;
    std::string        _name;
    CLHEP::HepRotation _rotation; // wrt to parent volume

  };

  std::ostream& operator<<(std::ostream& os, const Tube& t);
}

#endif/*GeomPrimitives_Tube_hh*/
