// Compute track parameter estimates by performing two separate
// linear regressions:  straight line in (ZX) projection, and
// "straight line with a kink" in (YZ) projection.
//
// Andrei Gaponenko, 2012

#include "ExtinctionMonitorFNAL/Reconstruction/inc/LinearRegression.hh"

#include <cmath>

#include "cetlib/exception.h"

#include "ExtinctionMonitorFNAL/Reconstruction/inc/Tracklet.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNAL.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALSensorStack.hh"

namespace mu2e {
  namespace ExtMonFNAL {

    //================================================================
    LinearRegression::LinearRegression(const ExtMon *extmon, double thetaScatterOnePlane)
      : extmon_(extmon)
      , scatterSigma2(extmon->nplanes())
      , zStart_(extmon->up().sensor_zoffset().back())
    {
      // We estimate pars at the end of the upstream stack. The scatter
      // there is zero by definition.  Each plane introduces a kink
      // with the given thetaScatterOnePlane RMS, which is multiplied by
      // the distance to the next plane to get scattering effect on position.
      // The effect is cumulative (quadratic sum).

      double current_z = extmon_->sensorCenterInExtMon(extmon_->nplanes()-1).z();
      scatterSigma2[extmon_->nplanes()-1] = 0;
      for(int i= extmon_->nplanes()-2; i>=0; --i) {
        const double new_z = extmon_->sensorCenterInExtMon(i).z();
        const double dz = new_z - current_z;
        current_z = new_z;
        const double segmentScatter2 = std::pow(thetaScatterOnePlane * dz, 2);
        scatterSigma2[i] = segmentScatter2 + scatterSigma2[i+1];
      }
    }


    //================================================================
    void LinearRegression::addCluster(LinearRegressionData *eqs, const ExtMonFNALRecoCluster& cl) {
      const double zKink = 0;
      const double D = zKink - zStart_;

      const CLHEP::Hep3Vector hitpos = extmon_->stackToExtMon_position(cl.position());
      const double di =  hitpos.z() - zStart_;

      const bool dn = (hitpos.z() < 0.);

      {
        const double cls2x = cl.xWidth() * extmon_->chip().xPitch() / sqrt(12.);
        const double sigma2x = cls2x + scatterSigma2[cl.plane()];
        eqs->zxA[0][0] += 1/sigma2x;
        eqs->zxA[0][1] += di/sigma2x;
        eqs->zxA[1][1] += di*di/sigma2x;

        eqs->zxRHS[0]  += hitpos.x()/sigma2x;
        eqs->zxRHS[1]  += hitpos.x()*di/sigma2x;
      }

      {
        const double cls2y = cl.yWidth() * extmon_->chip().yPitch() / sqrt(12.);
        const double sigma2y = cls2y + scatterSigma2[cl.plane()];

        eqs->yzA[0][0] += 1/sigma2y;
        eqs->yzA[0][1] += di/sigma2y;
        eqs->yzA[1][1] += di*di/sigma2y;

        if(dn) {
          eqs->yzA[0][2] += (di - D)/sigma2y;
          eqs->yzA[1][2] += di*(di - D)/sigma2y;
          eqs->yzA[2][2] += (di - D)*(di - D)/sigma2y;
        }

        eqs->yzRHS[0]  += hitpos.y()/sigma2y;
        eqs->yzRHS[1]  += hitpos.y()*di/sigma2y;
        if(dn) {
          eqs->yzRHS[2]  += hitpos.y()*(di - D)/sigma2y;
        }
      }
    }

    //================================================================
    ExtMonFNALTrkParam LinearRegression::estimatePars(const Tracklet& up, const Tracklet& dn) {
      LinearRegressionData eqs;

      addCluster(&eqs, *up.lastCluster);
      addCluster(&eqs, *up.firstCluster);
      for(unsigned i=0; i<up.middleClusters.size(); ++i) {
        addCluster(&eqs, *up.middleClusters[i]);
      }

      addCluster(&eqs, *dn.lastCluster);
      addCluster(&eqs, *dn.firstCluster);
      for(unsigned i=0; i<dn.middleClusters.size(); ++i) {
        addCluster(&eqs, *dn.middleClusters[i]);
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
