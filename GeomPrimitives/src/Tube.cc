//
// "Typical" Tube object
//
// $Id: Tube.cc,v 1.2 2012/04/05 18:43:07 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/05 18:43:07 $
//
// Original author KLG
//

#include <string>

#include "CLHEP/Vector/ThreeVector.h"

#include "GeomPrimitives/inc/Tube.hh"

namespace mu2e {

  Tube::Tube(double rIn, double rOut, double halfLength, double phi0, double phiMax,
             std::string const & materialName, CLHEP::Hep3Vector const & originInMu2e)
    :
    _params(rIn,
            rOut,
            halfLength,
            phi0,
            phiMax),
    _originInMu2e(originInMu2e),
    _materialName(materialName)
  {};

  Tube::Tube(std::string const & materialName, CLHEP::Hep3Vector const & originInMu2e,
       double rIn, double rOut, double halfLength, double phi0, double phiMax)
    :
    _params(rIn,
            rOut,
            halfLength,
            phi0,
            phiMax),
    _originInMu2e(originInMu2e),
    _materialName(materialName)
  {};

  const Tube Tube::UNINITIALIZED("", CLHEP::Hep3Vector(), 0, 0, 0);
  
}
