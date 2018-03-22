#ifndef Mu2eG4_constructStoppingTarget_hh
#define Mu2eG4_constructStoppingTarget_hh
//
// Free function to construct the stopping targets.
//
// $Id: constructStoppingTarget.hh,v 1.5 2011/08/04 18:52:25 genser Exp $
// $Author: genser $
// $Date: 2011/08/04 18:52:25 $
//
// Original author Peter Shanahan
//
// Notes:

#include "G4Helper/inc/VolumeInfo.hh"

namespace mu2e{

    class SimpleConfig;

    VolumeInfo constructStoppingTarget( VolumeInfo   const& mother,
                                      SimpleConfig const& config );
    
}  // end namespace mu2e

#endif /* Mu2eG4_constructStoppingTarget_hh */
