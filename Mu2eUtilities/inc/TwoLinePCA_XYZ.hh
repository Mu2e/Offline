#ifndef Mu2eUtilities_TwoLinePCA_XYZ_hh
#define Mu2eUtilities_TwoLinePCA_XYZ_hh
//
// Given two lines in 3D, compute the distance of closest
// approach between the two lines.  The lines are
// specified in point-slope form.
//
//
// Original author Rob Kutschke
//
// There is special treatment for lines that are almost parallel.
// Lines are defined as close to parallel if sin**2 of the angle
// between the two tracks is less than a given cut, which has a
// default value specified in the declaration of the constructor.
//
// Constructor arguments:
// p1 - a point on line 1.
// t1 - a vector in the direction of line 1.
// p2 - a point on line 2.
// tc - a vector in the direction of line 2.
//
// parallelcut - if the sin**2 of the angle between the lines
//               is less than this cut, then the lines
//               are declared parallel and the
//
//

#include "Offline/DataProducts/inc/GenVector.hh"

namespace mu2e {

  class TwoLinePCA_XYZ{

  public:
    TwoLinePCA_XYZ( XYZVectorF const& p1,
                XYZVectorF const& t1,
                XYZVectorF const& p2,
                XYZVectorF const& t2,
                double closeToParallelCut = 1.e-8
                );
    ~TwoLinePCA_XYZ();

    // Accessors for the 3d and 2d DCA.
    double dca()   const { return _dca; }
    double dca2d() const { return _dca2d; }
    int LRambig()   const { return _LRambig; }
    // Accessors for the endpoints on the line-segment of closest approach.
    XYZVectorF const& point1() const { return _pca1; }
    XYZVectorF const& point2() const { return _pca2; }

    // "Local z" of the pca.  ( Length along the straw from mid point to pca ).
    double s1() const { return _s1; }
    double s2() const { return _s2; }

    // Were the lines treated as parallel?
    bool closeToParallel() const { return _closeToParallel; }

  private:

    // A point on each track and a unit vector in the direction of each track.
    XYZVectorF _p1, _t1;
    XYZVectorF _p2, _t2;

    // Distances along each track from the points p1 or p2 to the points of closest approach.
    double _s1, _s2;

    // Points on each line at that terminate the line-segment of closest approach.
    XYZVectorF _pca1, _pca2;

    // Distance of closest approach in 3d and 2d(projection onto xy plane).
    double _dca;
    double _dca2d;

    // Did we trigger the close to parallel cut.
    bool _closeToParallel;

    // The value of that cut.
    double _cut;
    int _LRambig;
  };

} // namespace mu2e

#endif /* Mu2eUtilities_TwoLinePCA_XYZ_hh */
