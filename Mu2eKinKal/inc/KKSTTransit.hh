//
// KKSTTransit: estimate the stopping target (ST) material traversed by a fit trajectory, including its extrapolation.
// The trajectory is integrated in time through the hollow cylinder bounding the ST foils (ST_Front, ST_Back, ST_Inner, ST_Outer).
// At each time the probability that the true position is inside that volume is computed from the trajectory position covariance,
// projected onto the plane transverse to the velocity (along-track errors only re-label the time, so they integrate out), and
// then onto the radial and z directions (their correlation is ignored).  The time integral of speed times that probability is
// the expected path length in the volume.  Multiplying by the foil fill factor (total foil volume / bounding volume) estimates the
// expected path length through foil material.  The energy loss is estimated from the total material path length of each arc.
// Each contiguous sweep of the trajectory through the volume is recorded as a separate TrkArcTransit.
// The ST geometry (KinKalGeom) and material (KKMaterial) come from the GeometryService.
//
#ifndef Mu2eKinKal_KKSTTransit_hh
#define Mu2eKinKal_KKSTTransit_hh
#include "Offline/RecoDataProducts/inc/TrkMaterialTransit.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/KinKalGeom/inc/KinKalGeom.hh"
#include "Offline/KinKalGeom/inc/KKMaterial.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "fhiclcpp/types/Atom.h"
#include "cetlib_except/exception.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>

namespace mu2e {
  struct KKSTTransitConfig {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<double> stepSize { Name("StepSize"), Comment("Time integration step (ns)"), 0.01 };
    fhicl::Atom<double> pad { Name("Pad"), Comment("Trajectory points further than this from the bounding volume are skipped (mm)"), 100.0 };
    fhicl::Atom<double> minProb { Name("MinProb"), Comment("Minimum probability of being inside the bounding volume for a point to belong to an arc"), 1.0e-6 };
  };

  class KKSTTransit {
    public:
      using Config = KKSTTransitConfig;
      explicit KKSTTransit(Config const& config) : dt_(config.stepSize()), pad_(config.pad()), pmin_(config.minProb()) {
        if(dt_ <= 0.0) throw cet::exception("RECO") << "KKSTTransit: StepSize must be positive";
      }
      // (re)load the ST geometry and material; must be called when the geometry is available (ie beginRun)
      void beginRun();
      // ST material transit estimate of a trajectory
      template <class KTRAJ> TrkMaterialTransit STTransit(KinKal::ParticleTrajectory<KTRAJ> const& ptraj) const;
    private:
      double dt_, pad_, pmin_; // configuration
      double zlo_ = 0.0, zhi_ = 0.0, rin_ = 0.0, rout_ = 0.0; // bounding volume
      double fill_ = 0.0; // foil volume fill factor of the bounding volume
      SurfaceId szlo_, szhi_; // bounding disks at low and high z
      bool hasInner_ = false;
      std::shared_ptr<MatEnv::DetMaterial> mat_; // foil material
  };

  inline void KKSTTransit::beginRun() {
    GeomHandle<KinKalGeom> kkg_h;
    GeomHandle<KKMaterial> kkmat_h;
    auto const& st = *kkg_h->ST();
    auto const& outer = st.outer();
    zlo_ = (outer.center() - outer.axis()*outer.halfLength()).Z();
    zhi_ = (outer.center() + outer.axis()*outer.halfLength()).Z();
    if(zlo_ > zhi_) std::swap(zlo_,zhi_);
    rout_ = outer.radius();
    hasInner_ = st.hasInner();
    rin_ = hasInner_ ? st.inner().radius() : 0.0;
    bool frontlow = st.front().center().Z() < st.back().center().Z();
    szlo_ = frontlow ? SurfaceId(SurfaceIdEnum::ST_Front) : SurfaceId(SurfaceIdEnum::ST_Back);
    szhi_ = frontlow ? SurfaceId(SurfaceIdEnum::ST_Back) : SurfaceId(SurfaceIdEnum::ST_Front);
    double foilvol(0.0);
    for(size_t ifoil=0; ifoil < st.foils().size(); ++ifoil){
      auto const& foil = st.foil(ifoil);
      foilvol += (foil.outerRadius()*foil.outerRadius() - foil.innerRadius()*foil.innerRadius())*st.foilThickness(ifoil);
    }
    double boundvol = (rout_*rout_ - rin_*rin_)*(zhi_ - zlo_);
    if(boundvol <= 0.0) throw cet::exception("RECO") << "KKSTTransit: empty ST bounding volume";
    fill_ = foilvol/boundvol;
    mat_ = kkmat_h->STMaterial();
    if(!mat_) throw cet::exception("RECO") << "KKSTTransit: no ST material";
  }

  template <class KTRAJ> TrkMaterialTransit KKSTTransit::STTransit(KinKal::ParticleTrajectory<KTRAJ> const& ptraj) const {
    TrkMaterialTransit retval(SurfaceId(SurfaceIdEnum::ST_Foils));
    if(!mat_) throw cet::exception("RECO") << "KKSTTransit: beginRun not called";
    static const double isq2 = 1.0/std::sqrt(2.0);
    // probability that x is below hi (above lo) given a gaussian error sig
    auto below = [](double x, double hi, double sig) { return sig > 0.0 ? 0.5*std::erfc((x-hi)*isq2/sig) : (x < hi ? 1.0 : 0.0); };
    auto above = [](double x, double lo, double sig) { return sig > 0.0 ? 0.5*std::erfc((lo-x)*isq2/sig) : (x > lo ? 1.0 : 0.0); };
    auto inrange = [](double x, double lo, double hi, double sig) {
      if(sig <= 0.0) return (x > lo && x < hi) ? 1.0 : 0.0;
      return std::max(0.0,0.5*(std::erf((hi-x)*isq2/sig) - std::erf((lo-x)*isq2/sig)));
    };
    auto range = ptraj.range();
    TrkArcTransit arc;
    bool inarc(false);
    SurfaceId slimit; // limiting surface at the latest point in an arc
    double mass = ptraj.mass();
    // close an arc: compute its material and energy loss
    auto close = [&](){
      arc.matlen_ = fill_*arc.vollen_;
      if(arc.matlen_ > 0.0 && arc.mom_ > 0.0){
        arc.eloss_ = -mat_->energyLoss(arc.mom_,arc.matlen_,mass); // DetMaterial energy loss is negative
        arc.elossrms_ = mat_->energyLossRMS(arc.mom_,arc.matlen_,mass);
      }
      retval.arcs_.push_back(arc);
      arc = TrkArcTransit();
      inarc = false;
    };
    double time = range.begin() + 0.5*dt_;
    while(time < range.end()){
      auto pos = ptraj.position3(time);
      double speed = ptraj.velocity(time).R();
      double rho = pos.Rho();
      // distance outside the padded volume: the trajectory cannot get there faster than speed
      double dout = std::max({zlo_ - pad_ - pos.Z(), pos.Z() - zhi_ - pad_, rho - rout_ - pad_, 0.0});
      double prob(0.0);
      if(dout <= 0.0){
        auto pstate = ptraj.stateEstimate(time);
        auto const& cov = pstate.stateCovariance();
        double cphi = rho > 0.0 ? pos.X()/rho : 1.0;
        double sphi = rho > 0.0 ? pos.Y()/rho : 0.0;
        // project the radial and z directions onto the plane perpendicular to the velocity
        auto vdir = pstate.velocity().Unit();
        KinKal::VEC3 rdir(cphi,sphi,0.0), zdir(0.0,0.0,1.0);
        auto proj = [&cov](KinKal::VEC3 const& d){
          double dv[3] = {d.X(),d.Y(),d.Z()};
          double var(0.0);
          for(int i=0;i<3;++i) for(int j=0;j<3;++j) var += dv[i]*cov(i,j)*dv[j];
          return std::sqrt(std::max(0.0,var));
        };
        double sigr = proj(rdir - vdir*vdir.Dot(rdir));
        double sigz = proj(zdir - vdir*vdir.Dot(zdir));
        prob = inrange(pos.Z(),zlo_,zhi_,sigz)*inrange(rho,rin_,rout_,sigr);
        if(prob > pmin_){
          // the limiting bounding surface is the one with the smallest one-sided probability
          double plim = above(pos.Z(),zlo_,sigz); slimit = szlo_;
          double p = below(pos.Z(),zhi_,sigz);
          if(p < plim){ plim = p; slimit = szhi_; }
          p = below(rho,rout_,sigr);
          if(p < plim){ plim = p; slimit = SurfaceId(SurfaceIdEnum::ST_Outer); }
          if(hasInner_){
            p = above(rho,rin_,sigr);
            if(p < plim){ plim = p; slimit = SurfaceId(SurfaceIdEnum::ST_Inner); }
          }
          if(!inarc){
            inarc = true;
            arc.tentry_ = time - 0.5*dt_;
            arc.sentry_ = time - 0.5*dt_ > range.begin() ? slimit : SurfaceId(); // the trajectory starts inside
            arc.mom_ = ptraj.momentum(time);
          }
          arc.vollen_ += prob*speed*dt_;
          arc.texit_ = time + 0.5*dt_;
          arc.sexit_ = slimit;
        }
      }
      if(prob <= pmin_ && inarc) close();
      time += (dout > 0.0 && speed > 0.0) ? std::max(dt_,dout/speed) : dt_;
    }
    if(inarc){
      arc.sexit_ = SurfaceId(); // the trajectory ends inside
      close();
    }
    return retval;
  }
}
#endif
