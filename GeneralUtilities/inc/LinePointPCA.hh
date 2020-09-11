#ifndef Mu2eUtilities_LinePointPCA_hh
#define Mu2eUtilities_LinePointPCA_hh
//
// Given a line and an external point, find the point on the line that is
// closest to the external point.
//
//
// Original author Rob Kutschke
//

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class LinePointPCA{

  public:
    LinePointPCA( CLHEP::Hep3Vector const& p,
                  CLHEP::Hep3Vector const& t,
                  CLHEP::Hep3Vector const& q
                  );

    // Accept compiler generated versions of the
    // d'tor, copy cc'tor and assignment operator.


    // The 3d and 2d distances of closest approach.
    double dca()   const { return _dca; }
    double dca2d() const { return _dca2d; }

    // The point on the line that is closest to _q.
    CLHEP::Hep3Vector const& pca() const { return _pca; }

    // The unit vector in the direction from _pca to _q.
    CLHEP::Hep3Vector const& unit() const { return _unit; }

  private:

    // The line, in point slope form.
    CLHEP::Hep3Vector _p, _t;

    // The external point.
    CLHEP::Hep3Vector _q;

    // The point on the line that is closest to _q.
    CLHEP::Hep3Vector _pca;

    // Unit vector from pca to external refernce point.
    CLHEP::Hep3Vector _unit;

    // Distance along the line from p to the pca.
    double _s;

    // Distance of closest approach in 3d and 2d(projection onto xy plane).
    double _dca;
    double _dca2d;

  };

} // namespace mu2e

#endif /* Mu2eUtilities_LinePointPCA_hh */
