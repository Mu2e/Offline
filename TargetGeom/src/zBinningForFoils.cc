//
// Free function to compute a convenient binning for a histogram that will show 
// the z positions of the target foils.  If the spacing of foils is uniform, then all
// foil centers are guaranteed to be on bin centers. 
//
// See additional details in the comments in the header file.
//
// $Id: zBinningForFoils.cc,v 1.1 2011/05/15 21:14:09 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2011/05/15 21:14:09 $
//
// Original author Rob Kutschke
//

#include "FWCore/Utilities/interface/Exception.h"

#include "TargetGeom/inc/Target.hh"
#include "TargetGeom/inc/zBinningForFoils.hh"

namespace mu2e{

  Binning zBinningForFoils( Target const& target, int nBinsDZ ){

    int nfoils = target.nFoils();

    // TargetMaker.cc was supposed to make sure that this cannot happen.  Oh well.
    if ( nfoils < 1 ){
      throw cms::Exception("GEOM")
        << "zBinningForFoils: using a stopping target with no foils!\n";
    }

    // Special case: set bin size to the foil thickness instead of basing it
    // on the spacing between foils.
    if ( nfoils == 1 ){
      double dz = (nBinsDZ+0.5)*target.foil(0).halfThickness();
      double z0 = target.foil(0).center().z();
      int nbins = 2*nBinsDZ + 1;
      return Binning( nbins, z0-dz, z0+dz );
    }

    // Z positions of the centers of the first and last targets.
    TargetFoil const& f1 = target.foil(0);
    TargetFoil const& f2 = target.foil(nfoils-1);
    double z1 = f1.center().z();
    double z2 = f2.center().z();
  
    // Compute bin width (dz) as advertised for the general case.
    double dz = (z2-z1)/(nfoils-1)/nBinsDZ;

    // Compute the final answer.
    int nbins = (nfoils-1)*nBinsDZ + 1 + 2*nBinsDZ;
    double lo = z1-(nBinsDZ+0.5)*dz;
    double hi = z2+(nBinsDZ+0.5)*dz;

    return Binning( nbins, lo, hi);
  }

} // namespace mu2e
