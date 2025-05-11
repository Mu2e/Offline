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

namespace mu2e {
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
      // ElementXing interface
      void updateReference(PTRAJ const& ptraj) override;
      void updateState(MetaIterConfig const& config,bool first) override;
      Parameters params() const override;
      double time() const override { return inter_.time_; }
      double transitTime() const override;
      KTRAJ const& referenceTrajectory() const override { return *reftrajptr_; }
      std::vector<MaterialXing>const&  matXings() const override { return mxings_; }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // specific accessors
      auto const& intersection() const { return inter_; }
      auto const& material() const { return mat_; }
      auto const& surfaceId() const { return sid_; }
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
      double pathlen = thick_/(inter_.norm_.Dot(inter_.pdir_));
      mxings_.emplace_back(mat_,pathlen);
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
      double pathlen = thick_/(inter_.norm_.Dot(inter_.pdir_));
      mxings_.emplace_back(mat_,pathlen);
      fparams_ = this->parameterChange(varscale_);
    }
  }

  template <class KTRAJ,class SURF> Parameters KKShellXing<KTRAJ,SURF>::params() const {
    return fparams_;
  }

  template <class KTRAJ,class SURF> double KKShellXing<KTRAJ,SURF>::transitTime() const {
    double pathlen = thick_/(inter_.norm_.Dot(inter_.pdir_));
    return pathlen/reftrajptr_->speed();
  }

  template <class KTRAJ,class SURF> void KKShellXing<KTRAJ,SURF>::print(std::ostream& ost,int detail) const {
    ost <<"Shell Xing time " << this->time();
    ost << std::endl;
  }

}
#endif
