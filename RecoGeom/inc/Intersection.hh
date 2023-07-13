//
//  Calculate the intersection point of a Trajectory with a surface
//  This must be specialized for every case (every pair of trajectory and surface)
//  original author: David Brown (LBNL) 2023
//
#ifndef RecoGeom_Intersection_hh
#define RecoGeom_Intersection_hh
#include "Offline/RecoGeom/inc/InterData.hh"
#include "Offline/RecoGeom/inc/Cylinder.hh"
#include "Offline/RecoGeom/inc/Plane.hh"
#include "Offline/RecoGeom/inc/Annulus.hh"
#include "Offline/RecoGeom/inc/Rectangle.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/CentralHelix.hh"
#include "KinKal/Trajectory/KinematicLine.hh"

namespace mu2e {
  namespace RecoGeom {
    template <class TRAJ, class SURF> class Intersection {
      public:
        // constructor performs the intersction
        explicit Intersection(TRAJ const& ktraj, SURF const& surface,double timehint,double tolerance=1e-8);
        auto const& interData() const { return data_; }
        auto const& trajectory() const { return ktraj_; }
        auto const& surface() const { return surf_; }
      private:
        InterData data_; // result of this intersction
        TRAJ const& ktraj_; // trajectory
        SURF const& surf_; // surface
        double tol_; // tolerance (mm)
    };
    //
    // specializations for the different trajector and surface types
    //
    template<> Intersection<KinKal::LoopHelix,RecoGeom::Cylinder>::Intersection(KinKal::LoopHelix const& ktraj, RecoGeom::Cylinder const& surface, double timehint,double tolerance) : ktraj_(ktraj), surf_(surface), tol_(tolerance) {
      // temporary: do nothing
      std::cout << "LoopHelix and Cylinder status " << data_.flag_ << std::endl;
    }
    //
    template<> Intersection<KinKal::CentralHelix,RecoGeom::Cylinder>::Intersection(KinKal::CentralHelix const& ktraj, RecoGeom::Cylinder const& surface, double timehint,double tolerance) : ktraj_(ktraj), surf_(surface), tol_(tolerance) {
      // temporary: do nothing
      std::cout << "CentralHelix and Cylinder status " << data_.flag_ << std::endl;
    }
    //
    template<> Intersection<KinKal::KinematicLine,RecoGeom::Cylinder>::Intersection(KinKal::KinematicLine const& ktraj, RecoGeom::Cylinder const& cyl, double timehint,double tolerance) : ktraj_(ktraj), surf_(cyl), tol_(tolerance) {
      // start at the hint time.  Since this is a linear trajectory, the answer should be exact.
      double ttest = timehint;
      auto pos = ktraj_.position3(ttest);
      auto dir = ktraj_.direction(ttest);
      Ray ray(dir,pos);
      double dist;
      auto iflag = cyl.intersect(ray,dist,tol_);
      if(iflag.hasAnyProperty(IntersectFlag::onsurface)){
        data_.flag_ = iflag;
        data_.pos_ = ray.position(dist);
        data_.norm_ = cyl.normal(data_.pos_);
        data_.pdir_ = dir;
        // calculate the time
        data_.time_ = timehint + dist/ktraj_.speed(ttest);
      }
    }
  }
}

#endif
