// Figure out which simulation stage corresponds to a given SimParticle
// and provide access to its associated PhysicalVolumeInfo-s.
// The March 2021 update of PhysicalVolumeInfoMultiCollection makes
// this helper's function almost trivial; it is kept to avoid
// breakage of existing code that relies on the helper.
//
// Andrei Gaponenko, 2013, 2021

#ifndef Mu2eUtilities_PhysicalVolumeMultiHelper_hh
#define Mu2eUtilities_PhysicalVolumeMultiHelper_hh

#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"

namespace mu2e {
struct PhysicalVolumeInfo;

  class PhysicalVolumeMultiHelper {
  public:
    typedef PhysicalVolumeInfoMultiCollection::size_type size_type;

    PhysicalVolumeMultiHelper(const PhysicalVolumeInfoMultiCollection& coll);

    // the volumes
    const PhysicalVolumeInfo& startVolume(const SimParticle& p) const;
    const PhysicalVolumeInfo& endVolume(const SimParticle& p) const;

  private:
    const PhysicalVolumeInfoMultiCollection *pi_;
  };

} // end of namespace mu2e

#endif /* Mu2eUtilities_PhysicalVolumeMultiHelper_hh */
