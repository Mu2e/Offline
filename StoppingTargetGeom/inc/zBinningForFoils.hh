#ifndef StoppingTargetGeom_zBinningForFoils_hh
#define StoppingTargetGeom_zBinningForFoils_hh

//
// Free function to compute a convenient binning for a histogram that will show
// the z positions of the target foils.  If the spacing of foils is uniform, then all
// foil centers are guaranteed to be on bin centers.
//
// The argument controls the bin size: it specifies how many bins to make between
// the centers of adjacent foils.  The binning definition returned includes extra
// space at both ends, equal to the distance between foils; so for a detector with N
// foils, the histogram covers N+1 foil spacings.
//
// The case of a target with one foil is handed differently.
// This uses the DetectorSystem coordinates, not the Mu2e system.
//
// $Id: zBinningForFoils.hh,v 1.1 2013/05/31 20:04:27 gandr Exp $
// $Author: gandr $
// $Date: 2013/05/31 20:04:27 $
//
// Original author Rob Kutschke

#include "GeneralUtilities/inc/Binning.hh"

namespace mu2e {

  class StoppingTarget;

  Binning zBinningForFoils( StoppingTarget const& target, int nBinsDZ );

} // namespace mu2e

#endif /* StoppingTargetGeom_zBinningForFoils_hh */
