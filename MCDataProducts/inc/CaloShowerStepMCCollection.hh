#ifndef MCDataProducts_CaloShowerStepMCCollection_hh
#define MCDataProducts_CaloShowerStepMCCollection_hh

//
// Define a type for a collection of CaloShowerStepMC objects.
//
// Original author Bertrand Echenard
//

#include <vector>
#include "MCDataProducts/inc/CaloShowerStepMC.hh"

namespace mu2e {
   typedef std::vector<mu2e::CaloShowerStepMC> CaloShowerStepMCCollection;
}

#endif 
