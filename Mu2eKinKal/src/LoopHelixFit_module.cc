//
// KinKal fit module for LoopHelix
//
// Original author D. Brown (LBNL) 11/18/20
//
#include "KinKal/Trajectory/LoopHelix.hh"
using KTRAJ= KinKal::LoopHelix;
#include "Offline/Mu2eKinKal/inc/HelixFit_module.hh"
namespace mu2e {
  class LoopHelixFit : public HelixFit {
    public:
      explicit LoopHelixFit(const GlobalSettings& settings) :
        HelixFit(settings,TrkFitFlag::KKLoopHelix) {}
      KTRAJ makeSeedTraj(HelixSeed const& hseed) const override;
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
    Parameters kkpars(pars,seedcov_);
    //  construct the seed trajectory (infinite initial time range)
    return KTRAJ(kkpars, mass_, charge_, bnom, TimeRange(tmin,tmax));
  }
}
DEFINE_ART_MODULE(mu2e::LoopHelixFit);
