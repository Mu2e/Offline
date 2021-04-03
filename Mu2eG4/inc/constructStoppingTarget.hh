#ifndef Mu2eG4_constructStoppingTarget_hh
#define Mu2eG4_constructStoppingTarget_hh
//
// Free function to construct the stopping targets.
//
//
// Original author Peter Shanahan
//
// Notes:

#include "Mu2eG4Helper/inc/VolumeInfo.hh"

namespace mu2e{

    class SimpleConfig;

    VolumeInfo constructStoppingTarget( VolumeInfo   const& mother,
                                      SimpleConfig const& config );
    
}  // end namespace mu2e

#endif /* Mu2eG4_constructStoppingTarget_hh */
