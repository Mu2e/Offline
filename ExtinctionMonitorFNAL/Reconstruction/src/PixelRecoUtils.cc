#include "ExtinctionMonitorFNAL/Reconstruction/inc/PixelRecoUtils.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"

//================================================================
bool mu2e::ExtMonFNAL::perfectSingleParticleEvent(const ExtMonFNALRecoClusterCollection& coll, unsigned nExtMonPlanes) {
  bool accepted = true;
  for(unsigned plane = 0; plane < nExtMonPlanes; ++plane) {
    if(coll.clusters(plane).size() != 1) {
      accepted = false;
      break;
    }
  }

  return accepted;
}

//================================================================
