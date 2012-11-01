// Compute position uncertainty for a cluster assuming that track
// parameters are estimated at the last ustream sensor and including
// both cluster size and multiple scattering.
//
// Andrei Gaponenko, 2012

#ifndef ExtinctionMonitorFNAL_Reconstruction_ClusterOnTrackPrecisionTool_hh
#define ExtinctionMonitorFNAL_Reconstruction_ClusterOnTrackPrecisionTool_hh

#include <vector>

namespace mu2e {

  class ExtMonFNALRecoCluster;

  namespace ExtMonFNAL {

    class ExtMon;

    struct ClusterOnTrackPrecision {
      double sigma2x;
      double sigma2y;
      ClusterOnTrackPrecision(double s2x, double s2y) : sigma2x(s2x), sigma2y(s2y) {}
    };

    //================================================================
    class ClusterOnTrackPrecisionTool {
    public:

      explicit ClusterOnTrackPrecisionTool(const ExtMon& extmon, double thetaScatterOnePlane);

      ClusterOnTrackPrecision clusterPrecision(const ExtMonFNALRecoCluster& cl);

      // In mu2e software geometry is not accessible at module construction.
      // Provide default ctr so that modules don't have to deal with pointers to us.
      ClusterOnTrackPrecisionTool(): xPitch_(), yPitch_() {}

    private:
      // Pixel size
      double xPitch_;
      double yPitch_;
      std::vector<double> scatterSigma2_;
    };

  } // namespace ExtMonFNAL
} // namespace mu2e

#endif/*ExtinctionMonitorFNAL_Reconstruction_ClusterOnTrackPrecisionTool_hh*/
