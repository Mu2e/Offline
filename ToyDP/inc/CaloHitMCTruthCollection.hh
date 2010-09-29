#ifndef ToyDP_CaloHitMCTruthCollection_hh
#define ToyDP_CaloHitMCTruthCollection_hh

//
// Define a type for a collection of CaloHitMCTruth objects.
//
// Original author Ivan Logashenko
//

#include <vector>

#include "ToyDP/inc/CaloHitMCTruth.hh"

namespace mu2e {
   typedef std::vector<mu2e::CaloHitMCTruth> CaloHitMCTruthCollection;
}

#endif
