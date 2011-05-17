#ifndef TrackerGeom_TubsParams_hh
#define TrackerGeom_TubsParams_hh

//
// The parameters of a TUBS
//
// $Id: TubsParams.hh,v 1.3 2011/05/17 15:41:37 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/17 15:41:37 $
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


#endif /* TrackerGeom_TubsParams_hh */
