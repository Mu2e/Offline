#ifndef MCDataProducts_StrawDigiMCCollection_hh
#define MCDataProducts_StrawDigiMCCollection_hh
//
// Define a type for a collection of StrawDigiMC objects.
//
// Original author David Brown, LBNL 
//

#include <vector>
#include "MCDataProducts/inc/StrawDigiMC.hh"

namespace mu2e {
   typedef std::vector<mu2e::StrawDigiMC> StrawDigiMCCollection;
}
#endif /* MCDataProducts_StrawDigiMCCollection_hh */
