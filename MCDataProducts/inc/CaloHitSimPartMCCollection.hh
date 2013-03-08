#ifndef MCDataProducts_CaloHitSimPartMCCollection_hh
#define MCDataProducts_CaloHitSimPartMCCollection_hh

//
// Define a type for a collection of StrawHitMCTruth objects.
//
// Original author Bertrand Echenard
//

#include <vector>

#include "MCDataProducts/inc/CaloHitSimPartMC.hh"

namespace mu2e {
   typedef std::vector<mu2e::CaloHitSimPartMC> CaloHitSimPartMCCollection;
}

#endif /* MCDataProducts_CaloHitSimPartMCCollection_hh */
