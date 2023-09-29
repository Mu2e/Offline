//
// KinKal fit module for LoopHelix
//
// Original author D. Brown (LBNL) 11/18/20
//
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/General/Vectors.hh"
using KTRAJ= KinKal::LoopHelix;
#include "Offline/Mu2eKinKal/inc/HelixFit_module.hh"
namespace mu2e {
  using KinKal::DVEC;
  using KinKal::VEC3;
  class LoopHelixFit : public HelixFit {
    public:
      explicit LoopHelixFit(const Parameters& settings) :
        HelixFit(settings,TrkFitFlag::KKLoopHelix) {}
      // parameter-specific functions
      KTRAJ makeSeedTraj(HelixSeed const& hseed,TimeRange const& trange,VEC3 const& bnom, int charge) const override;
      bool goodFit(KKTRK const& ktrk) const override;
      virtual ~LoopHelixFit() {}
  };

  KTRAJ LoopHelixFit::makeSeedTraj(HelixSeed const& hseed,TimeRange const& trange,VEC3 const& bnom, int charge) const {
    auto const& helix = hseed.helix();
    DVEC pars;
    double psign = copysign(1.0,-charge*bnom.Z());
    pars[KTRAJ::rad_] = helix.radius()*psign;
    pars[KTRAJ::lam_] = helix.lambda()*fdir_.dzdt();
    pars[KTRAJ::cx_] = helix.centerx();
    pars[KTRAJ::cy_] = helix.centery();
    pars[KTRAJ::phi0_] = helix.fz0()+psign*M_PI_2;
    pars[KTRAJ::t0_] = hseed.t0().t0();
    // create the initial trajectory
    KinKal::Parameters kkpars(pars,seedcov_);
    //  construct the seed trajectory (infinite initial time range)
    return KTRAJ(kkpars, mass_, charge, bnom, trange);
  }

  bool LoopHelixFit::goodFit(KKTRK const& ktrk) const {
    // require physical consistency: fit can succeed but the result can have changed charge or helicity
    bool retval = ktrk.fitStatus().usable() &&
      ktrk.fitTraj().front().parameterSign()*ktrk.seedTraj().front().parameterSign() > 0 &&
      ktrk.fitTraj().front().helicity()*ktrk.seedTraj().front().helicity() > 0;
    // also check that the fit is inside the physical detector volume.  Test where the StrawHits are
    if(retval){
      for(auto const& shptr : ktrk.strawHits()) {
        if(shptr->active() && !Mu2eKinKal::inDetector(ktrk.fitTraj().position3(shptr->time()))){
          retval = false;
          break;
        }
      }
    }
    // test that the trajectory is inside the DS
    if(retval){
      static unsigned ntimes(100);
      double dt = ktrk.fitTraj().range().range()/(ntimes-1);
      for(unsigned it=0;it< ntimes; ++it) {
        double ttest = ktrk.fitTraj().range().begin() + it*dt;
        auto tpos = ktrk.fitTraj().position3(ttest);
        if(!ktrk.bfield().inRange(tpos)){
          retval = false;
          break;
        }
      }
    }
    return retval;
  }
}
DEFINE_ART_MODULE(mu2e::LoopHelixFit)
