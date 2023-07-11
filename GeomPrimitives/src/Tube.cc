//
// "Typical" Tube object
//
//
// Original author KLG
//

#include <string>

#include "Offline/GeomPrimitives/inc/Tube.hh"

namespace mu2e {

  Tube::Tube(double rIn, double rOut, double halfLength, double phi0, double phiMax,
             std::string const & materialName, CLHEP::Hep3Vector const & originInMu2e,
             CLHEP::HepRotation const & rotation )
    :
    _params(rIn,
            rOut,
            halfLength,
            phi0,
            phiMax),
    _originInMu2e(originInMu2e),
    _materialName(materialName),
    _rotation(rotation)
  {}

  Tube::Tube(double rIn, double rOut, double halfLength, CLHEP::Hep3Vector const & originInMu2e,
             CLHEP::HepRotation const & rotation,
             double phi0, double phiMax, std::string const & materialName )
    :
    _params(rIn,
            rOut,
            halfLength,
            phi0,
            phiMax),
    _originInMu2e(originInMu2e),
    _materialName(materialName),
    _rotation(rotation)
  {}

  Tube::Tube(std::string const & materialName, CLHEP::Hep3Vector const & originInMu2e,
             double rIn, double rOut, double halfLength, double phi0, double phiMax,
             CLHEP::HepRotation const & rotation )
    :
    _params(rIn,
            rOut,
            halfLength,
            phi0,
            phiMax),
    _originInMu2e(originInMu2e),
    _materialName(materialName),
    _rotation(rotation)
  {}

  std::ostream& operator<<(std::ostream& os, const Tube& t) {
    os<<"Tube: " << t.getTubsParams()
      <<", originInMu2e="<<t.originInMu2e()
      <<", phi0="<<t.phi0()
      <<", phiMax="<<t.phiMax()
      <<", rotation="<<t.rotation()
      <<", materialName="<<t.materialName();
    return os;
  }

}
