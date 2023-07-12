#include "Offline/RecoGeom/inc/Cylinder.hh"
#include "Math/VectorUtil.h"
using namespace ROOT::Math::VectorUtil;
namespace mu2e {
  namespace RecoGeom {
    bool Cylinder::onSurface(XYZVectorD const& point, double tol) const {
      auto rvec = point - center_;
      auto pvec = PerpVector(rvec,axis_);
      bool retval = fabs(pvec.R()-radius_) < tol;
      if(retval) retval = fabs(rvec.Dot(axis_)) - halflen_ < tol;
      return retval;
    }

    IntersectFlag Cylinder::intersect(Ray const& ray,double& dist, double tol) const {
      IntersectFlag retval;
      double ddot = ray.dir_.Dot(axis_);
      double alpha = (1.0 - ddot*ddot); // always positive
      // make sure the ray isn't co-linear
      if(alpha > tol/radius_){
        auto rvec = ray.start_ - center_;
        double sdot = rvec.Dot(axis_);
        double beta = sdot*ddot - rvec.Dot(ray.dir_);
        double gamma = rvec.Mag2() - sdot*sdot - radius2_;
        double beta2 = beta*beta;
        double ag = alpha*gamma;
        // make sure there's a solution
        if(beta2 > ag){
          double delta = sqrt(beta2 - ag);
          // choose smallest positive solution
          if(beta > delta){
            retval.merge(IntersectFlag::onsurface);
            dist = (beta - delta)/alpha;
          } else if( beta + delta > 0) {
            retval.merge(IntersectFlag::onsurface);
            dist = (beta + delta)/alpha;
          }
        }
      }
      // test that the point is on the cylinder
      if(retval.hasAllProperties(IntersectFlag::onsurface)){
        auto point = ray.position(dist);
        if(onSurface(point,tol))retval.merge(IntersectFlag::inbounds);
      }
      return retval;
    }
  }
}

std::ostream& operator <<(std::ostream& ost, mu2e::RecoGeom::Cylinder const& cyl) {
  ost << "Cylinder with center = " << cyl.center() << " , axis " << cyl.axis() << " radius " << cyl.radius() << " half-length " << cyl.halfLength();
  return ost;
}
