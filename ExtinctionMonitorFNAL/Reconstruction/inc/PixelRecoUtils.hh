#ifndef ExtinctionMonitorFNAL_Reconstruction_PixelRecoUtils_hh
#define ExtinctionMonitorFNAL_Reconstruction_PixelRecoUtils_hh

namespace mu2e {
  class ExtMonFNALRecoClusterCollection;

  namespace ExtMonFNAL {
    // Debugging aid
    bool perfectSingleParticleEvent(const ExtMonFNALRecoClusterCollection& coll, unsigned nExtMonPlanes);
  }

} // namespace mu2e

#endif/*ExtinctionMonitorFNAL_Reconstruction_PixelRecoUtils_hh*/
