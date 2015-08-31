// Compute track parameter estimates by performing two separate
// linear regressions:  straight line in (ZX) projection, and
// "straight line with a kink" in (YZ) projection.
//
// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Reconstruction/inc/ClusterOnTrackPrecisionTool.hh"

#include <cmath>

#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPlaneStack.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPixelChip.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    ClusterOnTrackPrecisionTool::ClusterOnTrackPrecisionTool(const ExtMon& extmon, double thetaScatterOnePlane)
      : xPitch_(extmon.chip().xPitch())
      , yPitch_(extmon.chip().yPitch())
      , scatterSigma2_(extmon.nplanes())
    {
      // We estimate pars at the end of the upstream stack. The scatter
      // there is zero by definition.  Each plane introduces a kink
      // with the given thetaScatterOnePlane RMS, which is multiplied by
      // the distance to the next plane to get scattering effect on position.
      // The effect is cumulative (quadratic sum).
      //
      // Could also add a displacement due to MS in air.  For that
      // it's better to add up rad lengths instead of adding up
      // displacement squares.

      double current_z = extmon.planeCenterInExtMon(extmon.nplanes()-1).z();
      scatterSigma2_[extmon.nplanes()-1] = 0;
      for(int i= extmon.nplanes()-2; i>=0; --i) {
        const double new_z = extmon.planeCenterInExtMon(i).z();
        const double dz = new_z - current_z;
        current_z = new_z;
        const double segmentScatter2 = std::pow(thetaScatterOnePlane * dz, 2);
        scatterSigma2_[i] = segmentScatter2 + scatterSigma2_[i+1];
      }
    }

    //================================================================
    ClusterOnTrackPrecision ClusterOnTrackPrecisionTool::clusterPrecision(const ExtMonFNALRecoCluster& cl) {

      const double cls2x = cl.xWidth() * xPitch_ / sqrt(12.);
      const double sigma2x = cls2x + scatterSigma2_[cl.plane()];

      const double cls2y = cl.yWidth() * yPitch_ / sqrt(12.);
      const double sigma2y = cls2y + scatterSigma2_[cl.plane()];

      return ClusterOnTrackPrecision(sigma2x, sigma2y);
    }

  } // namespace ExtMonFNAL
} // namespace mu2e
