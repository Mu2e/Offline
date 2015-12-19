#ifndef MCDataProducts_CaloDigiMCCollection_hh
#define MCDataProducts_CaloDigiMCCollection_hh

//
// Define a type for a collection of CaloDigiMC objects.
//
// Original author G. Pezzullo
//

#include <vector>

#include "MCDataProducts/inc/CaloDigiMC.hh"

namespace mu2e {
   typedef std::vector<mu2e::CaloDigiMC> CaloDigiMCCollection;
}

#endif /* MCDataProducts_CaloDigiMCCollection_hh */
