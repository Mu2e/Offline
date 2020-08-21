#ifndef AlgorithmID_AlgorithmIDCollection_hh
#define AlgorithmID_AlgorithmIDCollection_hh

//
// Define a type for a collection of AlgorithmID objects.
//
//
// Original author Vadim Rusu
//

#include <vector>

#include "RecoDataProducts/inc/AlgorithmID.hh"

namespace mu2e {
   typedef std::vector<mu2e::AlgorithmID> AlgorithmIDCollection;
}

#endif 
