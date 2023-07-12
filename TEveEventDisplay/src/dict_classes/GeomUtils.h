#ifndef GeomUtils_h
#define GeomUtils_h

//C++
#include <vector>
using namespace CLHEP;

namespace mu2e{
  inline constexpr double pointmmTocm(double mm){ return mm/10; }
  inline void hep3vectormmTocm(CLHEP::Hep3Vector &vector){vector.set(vector.x()/10, vector.y()/10, vector.z()/10);}

}
#endif
