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

namespace mu2e {
  template <class KTRAJ,class SURF> class KKShellXing : public ElementXing<KTRAJ> {
    public:
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using EXING = ElementXing<KTRAJ>;
      using PCA = PiecewiseClosestApproach<KTRAJ,SensorLine>;
      using CA = ClosestApproach<KTRAJ,SensorLine>;
      using SURFPTR = std::shared_ptr<SURF>;
      // construct from a surface, material, intersection, and transverse thickness
      KKShellXing(SURFPTR surface, MatEnv::DetMaterial const& mat, Intersection inter, double thickness)
      virtual ~KKShellXing() {}
      // ElementXing interface
      void updateReference(PTRAJ const& ptraj) override;
      void updateState(MetaIterConfig const& config,bool first) override;
      Parameters params() const override;
      double time() const override { return inter_.time(); }
      KTRAJ const& referenceTrajectory() const override { return *reftrajptr_; }
      std::vector<MaterialXing>const&  matXings() const override { return mxings_; }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // accessors
    private:
      SURFPTR surf_; // surface
      MatEnv::DetMaterial const& mat_;
      Intersection inter_; // most recent intersection
      KTRAJPTR reftrajptr_; // reference trajectory
      std::vector<MaterialXing> mxings_; // material xing
      double thick_; // shell thickness
      double tol_; // tolerance for intersection
      double varscale_; // variance scale, for annealing
      Parameters fparams_; // 1st-order parameter change for forwards time
  };

  template <class KTRAJ> KKShellXing<KTRAJ>::KKShellXing(SURFPTR surface, MatEnv::DetMaterial const& mat, Intersection inter,
      KTTRAJPTR reftrajptr, double thickness, double tol) :
    surf_(surface), mat_(mat), inter_(inter), reftrajptr_(reftrajptr), thick_(thick),tol_(tol),
    varscale_(1.0)
  {}

  template <class KTRAJ> void KKShellXing<KTRAJ>::updateReference(PTRAJ const& ptraj) {
    // re-intersect with the surface, taking the current time as start and range from the current piece (symmetrized)
    double delta = 0.5*reftraj->range().range();
    TimeRange irange(inter_.time_-delta, inter_.time_+delta);
    inter_ = intersect(ptraj, &surf_, irange,tol_);
    reftrajptr_ = ptraj.nearestTraj(inter_.time_);
  }

  template <class KTRAJ> void KKShellXing<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    // reset. This assumes we're off the surface
    fparams_ = Parameters();
    mxings_.clear();
    // check if we are on the surface
    if(inter_.onsurface_ && inter_.inbounds_){
      // compute the path length
      double pathlen = thick_/(inter_.norm_.Dot(inter_.pdir_);

    }
  }

  template <class KTRAJ> Parameters KKShellXing<KTRAJ>::params() const {
    return fparams_;
  }

  template <class KTRAJ> void KKShellXing<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost <<"Shell Xing time " << this->time();
    ost << std::endl;
  }

}
#endif
