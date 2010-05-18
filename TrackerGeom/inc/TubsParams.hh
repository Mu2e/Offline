#ifndef TUBSPARAMS_HH
#define TUBSPARAMS_HH

//
// The parameters of a TUBS
//
// $Id: TubsParams.hh,v 1.2 2010/05/18 20:29:06 kutschke Exp $
// $Author: kutschke $
// $Date: 2010/05/18 20:29:06 $
//
// Original author Rob Kutschke
//

#include <ostream>
#include <cmath>

#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e {

  struct TubsParams{
    double innerRadius;
    double outerRadius;
    double zHalfLength;
    double phi0;
    double phiMax;
    
    TubsParams( double innerRadius_,
                double outerRadius_,
                double zHalfLength_,
                double phi0_   = 0.,
                double phiMax_ = CLHEP::twopi):
      innerRadius(innerRadius_),
      outerRadius(outerRadius_),
      zHalfLength(zHalfLength_),
      phi0(phi0_),
      phiMax(phiMax_){
    }
      
    
  };

  inline std::ostream& operator<<(std::ostream& ost,
                                  const TubsParams& tp ){
    ost << "("
        << tp.innerRadius << " "
        << tp.outerRadius << " "
        << tp.zHalfLength << " "
        << tp.phi0        << " "
        << tp.phiMax      << ")"
        << " )";
    return ost;
  }

}


#endif
