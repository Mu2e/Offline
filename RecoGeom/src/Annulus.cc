#include "Offline/RecoGeom/inc/Annulus.hh"
#include "Math/VectorUtil.h"
using namespace ROOT::Math::VectorUtil;
namespace mu2e {
  namespace RecoGeom {
    bool Annulus::onSurface(XYZVectorD const& point, double tol) const {
      bool retval = Plane::onSurface(point,tol);
      if(retval){
        auto pvec = PerpVector(point - center(),normal());
        retval = onAnnulus(pvec.R());
      }
      return retval;
    }

    bool Annulus::intersect(Ray const& ray,double& dist, double tol) const {
      bool retval = Plane::intersect(ray,dist,tol);
      if(retval){
        retval = onSurface(ray.position(dist),tol);
      }
      return retval;
    }
  }
}
std::ostream& operator <<(std::ostream& ost, mu2e::RecoGeom::Annulus const& ann) {
  mu2e::RecoGeom::Plane const& plane = static_cast<mu2e::RecoGeom::Plane const&>(ann);
  ost << "Annulus in " << plane << " Inner, Outer radii " << ann.innerRadius() << " , " << ann.outerRadius();
  return ost;
}
