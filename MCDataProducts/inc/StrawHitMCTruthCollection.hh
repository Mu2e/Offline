#ifndef MCDataProducts_StrawHitMCTruthCollection_hh
#define MCDataProducts_StrawHitMCTruthCollection_hh

//
// Define a type for a collection of StrawHitMCTruth objects.
//
// Original author Ivan Logashenko
//

#include <vector>

#include "MCDataProducts/inc/StrawHitMCTruth.hh"

namespace mu2e {
   typedef std::vector<mu2e::StrawHitMCTruth> StrawHitMCTruthCollection;
}

#endif /* MCDataProducts_StrawHitMCTruthCollection_hh */
