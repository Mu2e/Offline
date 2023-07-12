#include "Offline/RecoGeom/inc/Plane.hh"
namespace mu2e {
  namespace RecoGeom {
    bool Plane::onSurface(XYZVectorD const& point, double tol) const {
      return fabs(norm_.Dot(point-center_)) < tol;
    }

    IntersectFlag Plane::intersect(Ray const& ray,double& dist, double tol) const {
      IntersectFlag retval;
      double ddir = norm_.Dot(ray.dir_);
      if(fabs(ddir)>tol){
        double pdist = norm_.Dot(center_ - ray.start_);
        dist = pdist/ddir;
        retval.merge(IntersectFlag::onsurface);
      }
      return retval;
    }
  }
}
std::ostream& operator <<(std::ostream& ost, mu2e::RecoGeom::Plane const& plane) {
  ost << "Plane with center " << plane.center() << " , normal " << plane.normal();
  return ost;
}
