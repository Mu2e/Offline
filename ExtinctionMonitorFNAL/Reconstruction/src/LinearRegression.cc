// Compute track parameter estimates by performing two separate
// linear regressions:  straight line in (ZX) projection, and
// "straight line with a kink" in (YZ) projection.
//
// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Reconstruction/inc/LinearRegression.hh"

#include <cmath>

#include "cetlib_except/exception.h"

#include "ExtinctionMonitorFNAL/Reconstruction/inc/Tracklet.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    LinearRegression::LinearRegression(const ExtMon *extmon, const ClusterOnTrackPrecisionTool& clTool)
      : extmon_(extmon)
      , clTool_(clTool)
      , zStart_(extmon->up().plane_zoffset().back())
    {}

    //================================================================
    void LinearRegression::addCluster(LinearRegressionData *eqs, const ExtMonFNALRecoCluster& cl) {
      const double zKink = 0;
      const double D = zKink - zStart_;

      const CLHEP::Hep3Vector hitpos = extmon_->stackToExtMon_position(cl.position());
      const double di =  hitpos.z() - zStart_;

      const ClusterOnTrackPrecision cp = clTool_.clusterPrecision(cl);

      const bool dn = (hitpos.z() < 0.);

      //----------------
      // ZX projection
      {
        eqs->zxA[0][0] += 1/cp.sigma2x;
        eqs->zxA[0][1] += di/cp.sigma2x;
        eqs->zxA[1][1] += di*di/cp.sigma2x;

        eqs->zxRHS[0]  += hitpos.x()/cp.sigma2x;
        eqs->zxRHS[1]  += hitpos.x()*di/cp.sigma2x;
      }

      //----------------
      // YZ projection
      {
        eqs->yzA[0][0] += 1/cp.sigma2y;
        eqs->yzA[0][1] += di/cp.sigma2y;
        eqs->yzA[1][1] += di*di/cp.sigma2y;

        if(dn) {
          eqs->yzA[0][2] += (di - D)/cp.sigma2y;
          eqs->yzA[1][2] += di*(di - D)/cp.sigma2y;
          eqs->yzA[2][2] += (di - D)*(di - D)/cp.sigma2y;
        }

        eqs->yzRHS[0]  += hitpos.y()/cp.sigma2y;
        eqs->yzRHS[1]  += hitpos.y()*di/cp.sigma2y;
        if(dn) {
          eqs->yzRHS[2]  += hitpos.y()*(di - D)/cp.sigma2y;
        }
      }
    }

    //================================================================
    ExtMonFNALTrkParam LinearRegression::estimatePars(const Tracklet& up, const Tracklet& dn) {
      LinearRegressionData eqs;

      addCluster(&eqs, *up.secondSeedCluster);
      addCluster(&eqs, *up.firstSeedCluster);
      for(unsigned i=0; i<up.addedClusters.size(); ++i) {
        addCluster(&eqs, *up.addedClusters[i]);
      }

      addCluster(&eqs, *dn.secondSeedCluster);
      addCluster(&eqs, *dn.firstSeedCluster);
      for(unsigned i=0; i<dn.addedClusters.size(); ++i) {
        addCluster(&eqs, *dn.addedClusters[i]);
      }

      // Solve the equations
      int ifail(0);

      eqs.zxA.invert(ifail);
      if(ifail) {
        throw cet::exception("DATA")<<"LinearRegression: Error inverting zxA matrix\n";
      }

      eqs.yzA.invert(ifail);
      if(ifail) {
        throw cet::exception("DATA")<<"LinearRegression: Error inverting yzA matrix\n";
      }

      const CLHEP::HepVector xpar(eqs.zxA * eqs.zxRHS);
      const CLHEP::HepVector ypar(eqs.yzA * eqs.yzRHS);

      // compute radius from the kink
      const double rinv = ypar[2]/(2.*extmon_->spectrometerMagnet().outerHalfSize()[2]);

      ExtMonFNALTrkParam res;
      res.setz0(zStart_);
      res.setposx(xpar[0]);
      res.setslopex(xpar[1]);
      res.setposy(ypar[0]);
      res.setslopey(ypar[1]);
      res.setrinv(rinv);
      return res;
    }

  } // namespace ExtMonFNAL
} // namespace mu2e
