//
// Given a line and an external point, find the point on the line that is
// closest to the external point.
//
// $Id: LinePointPCA.cc,v 1.1 2009/11/11 14:40:00 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/11 14:40:00 $
//
// Original author Rob Kutschke
//

#include "Mu2eUtilities/inc/LinePointPCA.hh"

using namespace std;

using CLHEP::Hep3Vector;

namespace mu2e{

  LinePointPCA::LinePointPCA( Hep3Vector const& p,
			      Hep3Vector const& t,
			      Hep3Vector const& q
			      ):
    _p(p), _q(q), _t(t){
    

    double s = _t.dot(_q-_p);

    _pca = _p + s*_t;
    
    Hep3Vector delta(_q-_pca);

    _unit  = delta.unit();
    _dca   = delta.mag();
    _dca2d = delta.perp();

  }

}
