// Compute track parameter estimates by performing two separate
// linear regressions:  straight line in (ZX) projection, and
// "straight line with a kink" in (YZ) projection.
//
// Andrei Gaponenko, 2012

#ifndef ExtinctionMonitorFNAL_Reconstruction_LinearRegression_hh
#define ExtinctionMonitorFNAL_Reconstruction_LinearRegression_hh

#include <vector>

#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"

#include "ExtinctionMonitorFNAL/Reconstruction/inc/ClusterOnTrackPrecisionTool.hh"

namespace mu2e {

  class ExtMonFNALRecoCluster;
  class ExtMonFNALTrkParam;

  namespace ExtMonFNAL {

    class Tracklet;
    class ExtMon;

    //================================================================
    class LinearRegression {
    public:
      ExtMonFNALTrkParam estimatePars(const Tracklet& up, const Tracklet& dn);

      // Ownership is not passed. The geometry object must be kept available while this
      // the regression object is in use.
      explicit LinearRegression(const ExtMon* extmon, const ClusterOnTrackPrecisionTool& clTool);

      // In mu2e software geometry is not accessible at module construction.
      // Provide default ctr so that modules don't have to deal with pointers to us.
      LinearRegression() : extmon_(), zStart_() {}

    private:

      // NB: We can't pre-compute inverted matrices because sigma
      // depends on the cluster size.
      const ExtMon *extmon_;
      ClusterOnTrackPrecisionTool clTool_;
      double zStart_;  // the coordinate where we want to measure track parameters

      struct LinearRegressionData {
        CLHEP::HepSymMatrix zxA;
        CLHEP::HepVector zxRHS;

        CLHEP::HepSymMatrix yzA;
        CLHEP::HepVector yzRHS;

        LinearRegressionData() : zxA(2), zxRHS(2), yzA(3), yzRHS(3) {}
      };

      void addCluster(LinearRegressionData *eqs, const ExtMonFNALRecoCluster& hit);

    };

  } // namespace ExtMonFNAL
} // namespace mu2e

#endif/*ExtinctionMonitorFNAL_Reconstruction_LinearRegression_hh*/
