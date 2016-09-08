#ifndef RecoDataProducts_KalSeedCollection_hh
#define RecoDataProducts_KalSeedCollection_hh

//
// Define a type for a collection of KalSeed objects.
// Original author David Brown
//

#include <vector>

#include "RecoDataProducts/inc/KalSeed.hh"

namespace mu2e {
  typedef std::vector<mu2e::KalSeed> KalSeedCollection;
}

#endif /* RecoDataProducts_KalSeedCollection_hh */
