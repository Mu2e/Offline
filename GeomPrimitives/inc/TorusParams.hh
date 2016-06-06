#ifndef GeomPrimitives_TorusParams_hh
#define GeomPrimitives_TorusParams_hh

//
// The parameters of a Torus, based on TubsParams by Rob Kutschke
//
//
// 
//

#include <cmath>
#include <ostream>

#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e {

  class TorusParams{

  public:

    TorusParams()
    {
      innRad = 0.;
      outRad = 0.;
      torRad = 0.;
      myPhi0 = 0.;
      myPMax = 0.;
    }

    TorusParams( double innerRadius,
                double outerRadius,
                double torusRadius,
                double phi0   = 0.,
                double phiMax = CLHEP::pi/2.0)
    {
      innRad = innerRadius;
      outRad = outerRadius;
      torRad = torusRadius;
      myPhi0 = phi0;
      myPMax = phiMax;
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    double innerRadius()  const { return innRad; }
    double outerRadius()  const { return outRad; }
    double torusRadius()  const { return torRad; }
    double phi0()         const { return myPhi0; }
    double phiMax()       const { return myPMax; }


  private:

    double innRad;
    double outRad;
    double torRad;
    double myPhi0;
    double myPMax;

  };

  inline std::ostream& operator<<(std::ostream& ost,
                                  const TorusParams& tp ){
    ost << "("
        << tp.innerRadius() << " "
        << tp.outerRadius() << " "
        << tp.torusRadius() << " "
        << tp.phi0()        << " "
        << tp.phiMax()      << ")"
        << " )";
    return ost;
  }

}


#endif /* GeomPrimitives_TorusParams_hh */
