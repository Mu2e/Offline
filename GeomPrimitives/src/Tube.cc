//
// "Typical" Tube object
//
// $Id: Tube.cc,v 1.3 2012/04/30 16:21:31 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/30 16:21:31 $
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

}
