// Each stage of multi-stage simulations may use a different geometry.
// This data product records the mapping between simulation stage
// (stored in each SimParticle) and the corresponding geometry info in
// PhysicalVolumeInfoSingleStage, so that  particle.simStage() is an index into the
// top level vector of PhysicalVolumeInfoMultiCollection.
//
// Andrei Gaponenko, 2013, 2021


#ifndef MCDataProducts_PhysicalVolumeInfoMultiCollection_hh
#define MCDataProducts_PhysicalVolumeInfoMultiCollection_hh

#include <vector>
#include <utility>
#include <ostream>

#include "cetlib/map_vector.h"

#include "MCDataProducts/inc/PhysicalVolumeInfo.hh"

namespace mu2e {
  typedef cet::map_vector<PhysicalVolumeInfo> PhysicalVolumeInfoSingleStage;
  typedef std::vector<PhysicalVolumeInfoSingleStage> PhysicalVolumeInfoMultiCollection;

  std::ostream& operator<<(std::ostream&, const PhysicalVolumeInfoSingleStage&);
  std::ostream& operator<<(std::ostream&, const PhysicalVolumeInfoMultiCollection&);
}

#endif /* MCDataProducts_PhysicalVolumeInfoMultiCollection_hh */
