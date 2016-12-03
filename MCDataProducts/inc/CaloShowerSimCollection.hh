#ifndef MCDataProducts_CaloShowerSimCollection_hh
#define MCDataProducts_CaloShowerSimCollection_hh

//
// Define a type for a collection of CaloShowerStep objects.
//
// Original author Bertrand Echenard
//

#include <vector>
#include "MCDataProducts/inc/CaloShowerSim.hh"

namespace mu2e {
   typedef std::vector<mu2e::CaloShowerSim> CaloShowerSimCollection;
}

#endif 
