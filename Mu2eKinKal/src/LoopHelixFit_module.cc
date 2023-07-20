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
      KTRAJ makeSeedTraj(HelixSeed const& hseed) const override;
      bool goodFit(KKTRK const& ktrk) const override;
      virtual ~LoopHelixFit() {}
  };

  KTRAJ LoopHelixFit::makeSeedTraj(HelixSeed const& hseed) const {
    // compute the magnetic field at the helix center.  We only want the z compontent, as the helix fit assumes B points along Z
    auto const& shelix = hseed.helix();
    double zmin = std::numeric_limits<float>::max();
    double zmax = std::numeric_limits<float>::min();
    double tmin = std::numeric_limits<float>::max();
    double tmax = std::numeric_limits<float>::min();
    auto const& hits = hseed.hits();
    for( auto const& hit : hits) {
      double zpos = hit.pos().z();
      zmin = std::min(zmin,zpos);
      zmax = std::max(zmax,zpos);
      tmin = std::min(tmin,(double)hit.correctedTime());
      tmax = std::max(tmax,(double)hit.correctedTime());
    }
    float zcent = 0.5*(zmin+zmax);
    VEC3 center(shelix.centerx(), shelix.centery(),zcent);
    auto bcent = kkbf_->fieldVect(center);
    VEC3 bnom(0.0,0.0,bcent.Z());
    // create a KTRAJ from the helix fit result, to seed the KinKal fit.  First, translate the parameters
    // Note the sign adjustments; RobustHelix is a purely geometric helix, with slightly different conventions
    DVEC pars;
    double psign = copysign(1.0,-charge_*bcent.Z());
    pars[KTRAJ::rad_] = shelix.radius()*psign;
    pars[KTRAJ::lam_] = shelix.lambda()*kkfit_.fitDirection().dzdt();
    pars[KTRAJ::cx_] = shelix.centerx();
    pars[KTRAJ::cy_] = shelix.centery();
    pars[KTRAJ::phi0_] = shelix.fz0()+psign*M_PI_2;
    pars[KTRAJ::t0_] = hseed.t0().t0();
    // create the initial trajectory
    KinKal::Parameters kkpars(pars,seedcov_);
    //  construct the seed trajectory (infinite initial time range)
    return KTRAJ(kkpars, mass_, charge_, bnom, TimeRange(tmin,tmax));
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
