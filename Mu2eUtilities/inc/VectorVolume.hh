#ifndef Mu2eUtilities_VectorVolume_hh
#define Mu2eUtilities_VectorVolume_hh
//
// Given a vector and a volume, gives the intersections between the two
//
// Original author Stefano Roberto Soleti
//
//
// Constructor arguments:

//
//

// CLHEP includes
#include "CLHEP/Vector/ThreeVector.h"
#include <vector>
#include "GeneralUtilities/inc/safeSqrt.hh"
using namespace std;

namespace mu2e {

  class VectorVolume{

  public:
    VectorVolume(CLHEP::Hep3Vector const& _orig,
                 CLHEP::Hep3Vector const& _dir,
                 float _xMin, float _xMax,
                 float _yMin, float _yMax,
                 float _zMin, float _zMax);
    ~VectorVolume();

    bool pointInBox(float x, float y,
                    float x0, float y0,
                    float x1, float y1);
    float distance(const CLHEP::Hep3Vector &u, const CLHEP::Hep3Vector &v);
    void calIntersections(std::vector<CLHEP::Hep3Vector> & intersections);

  private:
    CLHEP::Hep3Vector _orig;
    CLHEP::Hep3Vector _dir;
    float _xMin;
    float _xMax;
    float _yMin;
    float _yMax;
    float _zMin;
    float _zMax;

  };

} // namespace mu2e

#endif /* Mu2eUtilities_VectorVolume_hh */
