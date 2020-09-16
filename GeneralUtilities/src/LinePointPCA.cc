//
// Given a line and an external point, find the point on the line that is
// closest to the external point.
//
//
// Original author Rob Kutschke
//

#include "GeneralUtilities/inc/LinePointPCA.hh"

using namespace std;

using CLHEP::Hep3Vector;

namespace mu2e{

  LinePointPCA::LinePointPCA( CLHEP::Hep3Vector const& p,
                              CLHEP::Hep3Vector const& t,
                              CLHEP::Hep3Vector const& q
                              ):
    _p(p), _t(t), _q(q) {


    double s = _t.dot(_q-_p);

    _pca = _p + s*_t;

    CLHEP::Hep3Vector delta(_q-_pca);

    _unit  = delta.unit();
    _dca   = delta.mag();
    _dca2d = delta.perp();

  }

}
