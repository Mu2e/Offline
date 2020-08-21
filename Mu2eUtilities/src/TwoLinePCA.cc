//
// Given two lines in 3D, compute the distance of closest
// approach between the two lines.  The lines are
// specified in point-slope form.
//
//
// Original author Rob Kutschke
//

#include <iostream>

#include "CLHEP/Vector/ThreeVector.h"


#include "Mu2eUtilities/inc/TwoLinePCA.hh"

using namespace std;

namespace mu2e {


  TwoLinePCA::TwoLinePCA( CLHEP::Hep3Vector const& p1,
                          CLHEP::Hep3Vector const& t1,
                          CLHEP::Hep3Vector const& p2,
                          CLHEP::Hep3Vector const& t2,
                          double cut
                          ):
    _p1(p1),
    _t1(t1.unit()),
    _p2(p2),
    _t2(t2.unit()),
    _s1(0.),
    _s2(0.),
    _dca(0.),
    _dca2d(0.),
    _closeToParallel(false),
    _cut(cut){

    // Cosine of angle between the two line directions.
    double c(_t1.dot(_t2));

    // Sine-squared corresponding to c.
    double sinsq(1.-c*c);

    // Check for (almost) parallel lines.
    if ( sinsq < cut ){
      _closeToParallel = true;
      _pca1 = p1;
      _pca2 = p2;

      //Hep3Vector diff(_pca1-_pca2);
      //_dca   = diff.mag();
      //_dca2d = diff.perp();

      // return;

    }
    // Normal case: lines far from parallel.
    else {

      CLHEP::Hep3Vector delta(_p1-_p2);
      double dDotT1 = delta.dot(_t1);
      double dDotT2 = delta.dot(_t2);

      _s1 =  (dDotT2*c-dDotT1)/sinsq;
      _s2 = -(dDotT1*c-dDotT2)/sinsq;

      _pca1 = _p1 + _t1*_s1;
      _pca2 = _p2 + _t2*_s2;
    }

    CLHEP::Hep3Vector diff = (_pca1-_pca2);
    _dca   = diff.mag();
    _dca2d = diff.perp();

  }

  TwoLinePCA::~TwoLinePCA(){
  }

} // end namespace mu2e
