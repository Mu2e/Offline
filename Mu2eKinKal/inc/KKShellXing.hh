#ifndef Mu2eKinKal_KKShellXing_hh
#define Mu2eKinKal_KKShellXing_hh
//
//  Describe the effects of a kinematic trajectory crossing a thin shell of material defined by a surface
//  Used in the kinematic Kalman fit
//
#include "KinKal/Detector/ElementXing.hh"
#include "KinKal/Geometry/Surface.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "Offline/DataProducts/inc/SurfaceId.hh"
#include <algorithm>
#include <cmath>

namespace mu2e {
  // Ratio of the unrestricted Bethe ionization mean to KinKal's energyLoss (the restricted/moyal loss), >=1.
  // We correct the eloss for thick passive materials (ie. concrete blocks) by scaling the crossing path length.
  // Side-effect: the multiple-scattering variance scales too (~3-8% larger pointing sigma at high p),
  // which does not touch the |p| residual. f>=1 always (we never reduce the loss).
  inline double betheCorrectionFactor(MatEnv::DetMaterial const& mat, double mom, double pathlen, double mass) {
    if(mom <= 0.0 || pathlen == 0.0) return 1.0;
    static constexpr double me = 0.510998950;               // electron mass [MeV]
    double const beta2 = mom*mom/(mom*mom + mass*mass);
    double const bg2   = (mom*mom)/(mass*mass);             // (beta*gamma)^2
    double const gamma = std::sqrt(mom*mom + mass*mass)/mass;
    double const xi    = mat.eloss_xi(std::sqrt(beta2), pathlen); // (K/2)(Z/A)*rho*|path|/beta^2  [MeV]
    double const I     = mat.eexc();                        // mean excitation energy [MeV]
    double const Tmax  = 2.0*me*bg2/(1.0 + 2.0*gamma*me/mass + (me/mass)*(me/mass));
    double const delta = mat.densityCorrection(bg2);        // same density-effect KinKal uses
    // unrestricted Bethe MEAN loss, as a NEGATIVE energy change (DetMaterial::energyLoss is negative for loss)
    double const meanChange = -xi*( std::log(2.0*me*bg2/I) + std::log(Tmax/I) - 2.0*beta2 - delta );
    double const el = mat.energyLoss(mom, pathlen, mass);   // moyalmean, negative
    if(el >= 0.0) return 1.0;
    double const f = meanChange/el;                         // |mean|/|moyal|, both negative -> positive
    return f > 1.0 ? f : 1.0;                               // only ever increase the loss
  }

  template <class KTRAJ,class SURF> class KKShellXing : public KinKal::ElementXing<KTRAJ> {
    public:
      using PTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using EXING = KinKal::ElementXing<KTRAJ>;
      using PCA = KinKal::PiecewiseClosestApproach<KTRAJ,SensorLine>;
      using CA = KinKal::ClosestApproach<KTRAJ,SensorLine>;
      using SURFPTR = std::shared_ptr<SURF>;
      // construct from a surface, material, intersection, and transverse thickness
      KKShellXing(SURFPTR surface, SurfaceId const& sid, MatEnv::DetMaterial const& mat, KinKal::Intersection inter, KTRAJPTR reftraj, double thickness, double tol);
      virtual ~KKShellXing() {}
      // clone op for reinstantiation
      KKShellXing(KKShellXing const& rhs) = default;
      std::shared_ptr< KinKal::ElementXing<KTRAJ> > clone(CloneContext& context) const override{
        auto rv = std::make_shared< KKShellXing<KTRAJ,SURF> >(*this);
        //auto ptr = context.get(reftrajptr_);
        auto ptr = std::make_shared<KTRAJ>(*reftrajptr_);
        rv->setReferenceTrajectoryPtr(ptr);
        return rv;
      };
      // ElementXing interface
      void updateReference(PTRAJ const& ptraj) override;
      void updateState(MetaIterConfig const& config,bool first) override;
      Parameters params() const override;
      double time() const override { return inter_.time_; }
      double transitTime() const override;
      KTRAJ const& referenceTrajectory() const override { return *reftrajptr_; }
      std::vector<MaterialXing>const&  matXings() const override { return mxings_; }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      bool active() const override { return mxings_.size() > 0; }
      // specific accessors
      auto const& intersection() const { return inter_; }
      auto const& material() const { return mat_; }
      auto const& surfaceId() const { return sid_; }
      // other accessors
      void setReferenceTrajectoryPtr(KTRAJPTR ptr){ reftrajptr_ = ptr; }
    private:
      SURFPTR surf_; // surface
      SurfaceId sid_; // surface Id
      MatEnv::DetMaterial const& mat_;
      KinKal::Intersection inter_; // most recent intersection
      KTRAJPTR reftrajptr_; // reference trajectory
      std::vector<KinKal::MaterialXing> mxings_; // material xing
      double thick_; // shell thickness
      double tol_; // tolerance for intersection
      double varscale_; // variance scale, for annealing
      KinKal::Parameters fparams_; // 1st-order parameter change for forwards time
  };

  template <class KTRAJ,class SURF> KKShellXing<KTRAJ,SURF>::KKShellXing(SURFPTR surface, SurfaceId const& sid, MatEnv::DetMaterial const& mat, KinKal::Intersection inter,
      std::shared_ptr<KTRAJ> reftrajptr, double thickness, double tol) :
    surf_(surface), sid_(sid), mat_(mat), inter_(inter), reftrajptr_(reftrajptr), thick_(thickness),tol_(tol),
    varscale_(1.0)
  {
    if(inter_.good()){
      // compute the path length
      double dotprod = std::max(1e-6,std::fabs(inter_.norm_.Dot(inter_.pdir_)));
      double pathlen = thick_/dotprod;
      // Scale the crossed path length so the (restricted/moyal) E loss becomes the unrestricted Bethe mean.
      double const f = betheCorrectionFactor(mat_, reftrajptr_->momentum(), pathlen, reftrajptr_->mass());
      mxings_.emplace_back(mat_, pathlen*f);
    }
  }

  template <class KTRAJ,class SURF> void KKShellXing<KTRAJ,SURF>::updateReference(KinKal::ParticleTrajectory<KTRAJ> const& ptraj) {
    auto const& neartrajptr = ptraj.nearestTraj(inter_.time_);
    if(neartrajptr != reftrajptr_){
      // traj has changed; try intersecting with this piece
      reftrajptr_ = neartrajptr;
      inter_ = KinKal::intersect(*reftrajptr_,*surf_, reftrajptr_->range(),tol_);
      // I should test if this is in range, but there's a missing piece in the ptraj intersection here that I don't know what direction,
      // to go, and there can be multiple intersections. Not clear if there's a good answer here TODO
    }
  }

  template <class KTRAJ,class SURF> void KKShellXing<KTRAJ,SURF>::updateState(KinKal::MetaIterConfig const& miconfig,bool first) {
    // reset. This assumes we're off the surface
    fparams_ = KinKal::Parameters();
    mxings_.clear();
    // check if we are on the surface; if so, create the xing
    if(inter_.good()){
      // compute the path length
      double dotprod = std::max(1e-6,std::fabs(inter_.norm_.Dot(inter_.pdir_)));
      double pathlen = thick_/dotprod;
      // same Bethe-mean path scale as the ctor (keeps the fit-side params consistent with the extrapolation)
      double const f = betheCorrectionFactor(mat_, reftrajptr_->momentum(), pathlen, reftrajptr_->mass());
      mxings_.emplace_back(mat_, pathlen*f);
      fparams_ = this->parameterChange(varscale_);
    }
  }

  template <class KTRAJ,class SURF> Parameters KKShellXing<KTRAJ,SURF>::params() const {
    return fparams_;
  }

  template <class KTRAJ,class SURF> double KKShellXing<KTRAJ,SURF>::transitTime() const {
    double dotprod = std::max(1e-6,fabs(inter_.norm_.Dot(inter_.pdir_)));
    double pathlen = thick_/dotprod;
    return pathlen/reftrajptr_->speed();
  }

  template <class KTRAJ,class SURF> void KKShellXing<KTRAJ,SURF>::print(std::ostream& ost,int detail) const {
    ost <<"Shell Xing time " << this->time();
    ost << std::endl;
  }

}
#endif
