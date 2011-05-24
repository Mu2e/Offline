#ifndef MCDataProducts_G4BeamlineInfoCollection_hh
#define MCDataProducts_G4BeamlineInfoCollection_hh

//
// Define a type for a collection of G4BeamlineInfo objects.
//

#include <vector>

#include "MCDataProducts/inc/G4BeamlineInfo.hh"

namespace mu2e {
   typedef std::vector<mu2e::G4BeamlineInfo> G4BeamlineInfoCollection;
}

#endif /* MCDataProducts_G4BeamlineInfoCollection_hh */
