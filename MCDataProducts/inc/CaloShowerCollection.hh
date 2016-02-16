#ifndef MCDataProducts_CaloShowerCollection_hh
#define MCDataProducts_CaloShowerCollection_hh

//
// Define a type for a collection of CaloShowerStepMC objects.
//
// Original author Bertrand Echenard
//

#include <vector>
#include "MCDataProducts/inc/CaloShower.hh"

namespace mu2e {
   typedef std::vector<mu2e::CaloShower> CaloShowerCollection;
}

#endif 
