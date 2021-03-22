// Andrei Gaponenko, 2013

#include <memory>

#include "cetlib/map_vector.h"
#include "cetlib_except/exception.h"

#include "Mu2eUtilities/inc/PhysicalVolumeMultiHelper.hh"

namespace mu2e {
struct PhysicalVolumeInfo;

  PhysicalVolumeMultiHelper::PhysicalVolumeMultiHelper(const PhysicalVolumeInfoMultiCollection& coll)
    : pi_(&coll)
  {}

  const PhysicalVolumeInfo&
  PhysicalVolumeMultiHelper::startVolume(const SimParticle& p) const {
    return (*pi_)[p.simStage()][cet::map_vector_key(p.startVolumeIndex())];
  }

  const PhysicalVolumeInfo&
  PhysicalVolumeMultiHelper::endVolume(const SimParticle& p) const {
    return (*pi_)[p.simStage()][cet::map_vector_key(p.endVolumeIndex())];
  }
}
