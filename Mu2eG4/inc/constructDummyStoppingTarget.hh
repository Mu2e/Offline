#ifndef Mu2eG4_constructDummyStoppingTarget_hh
#define Mu2eG4_constructDummyStoppingTarget_hh
//
// Free function to construct a placeholder for the stopping target.
// Useful for some low detail graphics.
//
// $Id: constructDummyStoppingTarget.hh,v 1.4 2011/08/04 18:52:25 genser Exp $
// $Author: genser $
// $Date: 2011/08/04 18:52:25 $
//
// Original author Rob Kutschke
//

#include "G4Helper/inc/VolumeInfo.hh"

namespace mu2e{

  class SimpleConfig;

  VolumeInfo constructDummyStoppingTarget( VolumeInfo   const& mother,
                                           SimpleConfig const& config );

}  // end namespace mu2e

#endif /* Mu2eG4_constructDummyStoppingTarget_hh */
