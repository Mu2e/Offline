#ifndef MCDataProducts_CaloShowerStepCollection_hh
#define MCDataProducts_CaloShowerStepCollection_hh

//
// Define a type for a collection of CaloShowerStep objects.
//
// Original author Bertrand Echenard
//

#include <vector>
#include "MCDataProducts/inc/CaloShowerStep.hh"

namespace mu2e {
   typedef std::vector<mu2e::CaloShowerStep> CaloShowerStepCollection;
}

#endif 
