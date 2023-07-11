#include "Offline/RecoGeom/inc/Rectangle.hh"
#include "Math/VectorUtil.h"
using namespace ROOT::Math::VectorUtil;
namespace mu2e {
  namespace RecoGeom {
    bool Rectangle::onSurface(XYZVectorD const& point, double tol) const {
      bool retval = Plane::onSurface(point,tol);
      if(retval){
        auto pvec = PerpVector(point - center(),normal());
        double udist = pvec.Dot(udir_);
        double vdist = pvec.Dot(vdir_);
        retval = onRectangle(udist,vdist);
      }
      return retval;
    }

    bool Rectangle::intersect(Ray const& ray,double& dist, double tol) const {
      bool retval = Plane::intersect(ray,dist,tol);
      if(retval){
        retval = onSurface(ray.position(dist),tol);
      }
      return retval;
    }
  }
}
std::ostream& operator <<(std::ostream& ost, mu2e::RecoGeom::Rectangle const& rect) {
  mu2e::RecoGeom::Plane const& plane = static_cast<mu2e::RecoGeom::Plane const&>(rect);
  ost << "Rectangle in " << plane << " U direction " << rect.uDirection() << " V direction " << rect.vDirection()
    << " U, V half-lengths " << rect.uHalfLength() << " , " << rect.vHalfLength();
  return ost;
}
