// Each stage of multi-stage simulations, in general, need its own
// PhysicalVolumeInfoCollection.   This data product records the mapping
// between simulation stages (identified by the sim particle number offset
// in the first element of the pair) and their PhysicalVolumeInfoCollections.
//
// Andrei Gaponenko, 2013


#ifndef MCDataProducts_PhysicalVolumeInfoMultiCollection_hh
#define MCDataProducts_PhysicalVolumeInfoMultiCollection_hh

#include <vector>
#include <utility>

#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"

namespace mu2e {
  typedef std::vector<std::pair<unsigned int, PhysicalVolumeInfoCollection> > PhysicalVolumeInfoMultiCollection;
}

#endif /* MCDataProducts_PhysicalVolumeInfoMultiCollection_hh */
