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
    template <class KTRAJ, class SURF> struct Intersection : public InterData {
      Intersection(KTRAJ const& ktraj, SURF const& surf,double tol) : ktraj_(ktraj), surf_(surf), tol_(tol) {}
      KTRAJ const& ktraj_; // trajectory of this intersection
      SURF const& surf_; // surface of this intersection
      double tol_; // tolerance used in this intersection
    };

    template <class KTRAJ, class SURF> struct Intersect {
      // perform the intersction.  Tolerance is in mm
      Intersection<KTRAJ, SURF> intersect( KTRAJ const& ktraj, SURF const& surface, double starttime ,double tolerance=1e-8) {
        return Intersection<KTRAJ, SURF>(ktraj,surface,tolerance);
      }
    };
    //
    // specializations for the different trajector and surface types
    //
    template<> struct Intersect<KinKal::LoopHelix,RecoGeom::Cylinder> {
      Intersection<KinKal::LoopHelix, RecoGeom::Cylinder> intersect( KinKal::LoopHelix const& ktraj, RecoGeom::Cylinder const& surface, double starttime ,double tolerance=1e-8) {
        Intersection<KinKal::LoopHelix, RecoGeom::Cylinder> retval(ktraj,surface,tolerance);
        // temporary: do nothing
        std::cout << "LoopHelix and Cylinder intersection NOT YET IMPLEMENTED" << std::endl;
        return retval;
      }
    };

    //
    // Line trajectory can provide an exact answer for generic surfaces
    //
    template<class SURF> struct Intersect<KinKal::KinematicLine,SURF> {
      Intersection<KinKal::KinematicLine,SURF> intersect(KinKal::KinematicLine const& ktraj, SURF const& surf, double starttime,double tolerance) {
        Intersection<KinKal::KinematicLine,SURF> retval(ktraj,surf,tolerance);
        auto pos = ktraj.position3(starttime);
        auto dir = ktraj.direction(starttime);
        Ray ray(dir,pos);
        double dist;
        retval.flag_ = surf.intersect(ray,dist,tolerance);
        if(retval.flag_.onsurface_){
          retval.pos_ = ray.position(dist);
          retval.norm_ = surf.normal(retval.pos_);
          retval.pdir_ = dir;
          // calculate the time
          retval.time_ = starttime + dist/ktraj.speed(starttime);
        }
        return retval;
      }
    };
  }
}

#endif
