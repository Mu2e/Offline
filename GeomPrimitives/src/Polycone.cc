//
// The parameters of a Polycone
//
// $Id: Polycone.cc,v 1.2 2012/03/29 19:07:00 gandr Exp $
// $Author: gandr $
// $Date: 2012/03/29 19:07:00 $
//
// Original author KLG
//

#include <cmath>
#include <string>
#include <algorithm>
#include <iterator>
#include <ostream>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "GeomPrimitives/inc/Polycone.hh"


namespace mu2e {

  Polycone::Polycone( unsigned int numZPlanes,
                      double zPlanes[],
                      double rInner[],
                      double rOuter[],
                      CLHEP::Hep3Vector const & originInMu2e,
                      std::string const & materialName,
                      double phiStart,
                      double phiTotal):
    _numZPlanes(numZPlanes),
    _zPlanes(),
    _rInner(),
    _rOuter(),
    _originInMu2e(originInMu2e),
    _materialName(materialName),
    _phiStart(phiStart),
    _phiTotal(phiTotal)
  {
    _zPlanes.reserve(_numZPlanes);
    std::copy (&zPlanes[0],&zPlanes[_numZPlanes],
               back_inserter(_zPlanes));
    _rInner.reserve(_numZPlanes);
    std::copy (&rInner[0],&rInner[_numZPlanes],
               back_inserter(_rInner));
    _rOuter.reserve(_numZPlanes);
    std::copy (&rOuter[0],&rOuter[_numZPlanes],
               back_inserter(_rOuter));
  };

  Polycone::Polycone( double phiStart,
                      double phiTotal,
                      unsigned int numZPlanes,
                      double zPlanes[],
                      double rInner[],
                      double rOuter[],
                      CLHEP::Hep3Vector const & originInMu2e,
                      std::string const & materialName) :
    _numZPlanes(numZPlanes),
    _zPlanes(),
    _rInner(),
    _rOuter(),
    _originInMu2e(originInMu2e),
    _materialName(materialName),
    _phiStart(phiStart),
    _phiTotal(phiTotal)
  {
    // FIXME one should factorize this
    _zPlanes.reserve(_numZPlanes);
    std::copy (&zPlanes[0],&zPlanes[_numZPlanes],
               back_inserter(_zPlanes));
    _rInner.reserve(_numZPlanes);
    std::copy (&rInner[0],&rInner[_numZPlanes],
               back_inserter(_rInner));
    _rOuter.reserve(_numZPlanes);
    std::copy (&rOuter[0],&rOuter[_numZPlanes],
               back_inserter(_rOuter));
  };

}
