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

  PhysicalVolumeMultiHelper::size_type
  PhysicalVolumeMultiHelper::iSimStage(SimParticleCollection::key_type key) const {
    size_type res = -1u;

    // There are only a few simulation stages.  Linear search is perhaps the most efficient way.
    for(int i = pi_->size() - 1; i >=0 ; --i) {
      if((*pi_)[i].first < key.asUint()) {
        res = i;
        break;
      }
    }

    if(res == -1u) {
      throw cet::exception("BADINPUTS")
        <<"PhysicalVolumeMultiHelper: simulation stage for key="<<key
        <<" is not in the collection\n";
    }

    return res;
  }

  const PhysicalVolumeInfo&
  PhysicalVolumeMultiHelper::startVolume(const SimParticle& p) const {
    return (*pi_)[iSimStage(p)].second[cet::map_vector_key(p.startVolumeIndex())];
  }

  const PhysicalVolumeInfo&
  PhysicalVolumeMultiHelper::endVolume(const SimParticle& p) const {
    return (*pi_)[iSimStage(p)].second[cet::map_vector_key(p.endVolumeIndex())];
  }


}
