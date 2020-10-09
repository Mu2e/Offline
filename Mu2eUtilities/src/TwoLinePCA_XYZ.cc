//
// Given two lines in 3D, compute the distance of closest
// approach between the two lines.  The lines are
// specified in point-slope form.
//
//
// Original author Rob Kutschke
//

#include <iostream>

#include "DataProducts/inc/XYZVec.hh"
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/DisplacementVector3D.h"

#include "Mu2eUtilities/inc/TwoLinePCA_XYZ.hh"

using namespace std;

namespace mu2e {


  TwoLinePCA_XYZ::TwoLinePCA_XYZ( XYZVec const& p1,
                          XYZVec const& t1,
                          XYZVec const& p2,
                          XYZVec const& t2,
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
    _cut(cut),
    _LRambig(0.)
     {

    // Cosine of angle between the two line directions.
    double c(_t1.Dot(_t2));

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

      XYZVec delta(_p1-_p2);
      double dDotT1 = delta.Dot(_t1);
      double dDotT2 = delta.Dot(_t2);

      _s1 =  (dDotT2*c-dDotT1)/sinsq;
      _s2 = -(dDotT1*c-dDotT2)/sinsq;

      _pca1 = _p1 + _t1*_s1;
      _pca2 = _p2 + _t2*_s2;
      _LRambig = _s2 > 0 ? 1 : -1;//int ambig_sign= ambig > 0 ? 1 : -1;
    }

    XYZVec diff = (_pca1-_pca2);
    _dca   = sqrt(diff.Mag2());
    _dca2d = sqrt(diff.Perp2());

  }

  TwoLinePCA_XYZ::~TwoLinePCA_XYZ(){
  }

} // end namespace mu2e
