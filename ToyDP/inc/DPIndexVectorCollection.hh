#ifndef ToyDP_DPIndexVectorCollection_hh
#define ToyDP_DPIndexVectorCollection_hh

//
// Define a type for a collection of DPIndexVector objects.
//
// Original author Ivan Logashenko
//

#include <vector>

#include "ToyDP/inc/DPIndexVector.hh"

namespace mu2e {
   typedef std::vector<mu2e::DPIndexVector> DPIndexVectorCollection;
}

#endif
