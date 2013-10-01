// Each stage of multi-stage simulations, in general, need its own
// PhysicalVolumeInfoCollection.   This data product records the mapping
// between simulation stages (identified by the sim particle number offset
// in the first element of the pair) and their PhysicalVolumeInfoCollections.
//
// Sim particle offsets (the first member of the pair) must be
// monotonically increasing with index in the collection.
//
// Andrei Gaponenko, 2013


#ifndef MCDataProducts_PhysicalVolumeInfoMultiCollection_hh
#define MCDataProducts_PhysicalVolumeInfoMultiCollection_hh

#include <vector>
#include <utility>

#include "cetlib/map_vector.h"

#include "MCDataProducts/inc/PhysicalVolumeInfo.hh"

namespace mu2e {
  typedef cet::map_vector<PhysicalVolumeInfo> PhysicalVolumeInfoSingleStage;
  typedef std::vector<std::pair<unsigned int, PhysicalVolumeInfoSingleStage> > PhysicalVolumeInfoMultiCollection;
}

#endif /* MCDataProducts_PhysicalVolumeInfoMultiCollection_hh */
