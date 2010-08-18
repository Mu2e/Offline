#ifndef ToyDP_StrawHitMCPtrCollection_hh
#define ToyDP_StrawHitMCPtrCollection_hh

//
// Define a type for a collection of StrawHitMCPtr objects.
//
// Original author Ivan Logashenko
//

#include <vector>

#include "ToyDP/inc/StrawHitMCPtr.hh"

namespace mu2e {
   typedef std::vector<mu2e::StrawHitMCPtr> StrawHitMCPtrCollection;
}

#endif
