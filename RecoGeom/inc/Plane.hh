//
//  Description of an unbounded plane in space.
//  original author: David Brown (LBN) 2023
//
#ifndef RecoGeom_Plane_hh
#define RecoGeom_Plane_hh
#include "Offline/RecoGeom/inc/Surface.hh"
namespace mu2e {
  namespace RecoGeom {
    class Plane : public Surface {
      public:
        // construct from necessary parameters
        Plane(XYZVectorD const& norm, XYZVectorD const& center) : norm_(norm.Unit()), center_(center){}
        // surface interface
        bool onSurface(XYZVectorD const& point, double tol=1e-8) const override;
        bool inBounds(XYZVectorD const& point, double tol=1e-8) const override { return true; }
        IntersectFlag intersect(Ray const& ray,double& dist, double tol=1e-8) const override;
        XYZVectorD normal(XYZVectorD const& point) const override { return norm_; }
        auto const& normal() const { return norm_; }
        auto const& center() const { return center_; }
      private:
        XYZVectorD norm_; // normal to the plane
        XYZVectorD center_; // point on the plane
    };
  }
}
std::ostream& operator <<(std::ostream& ost, mu2e::RecoGeom::Plane const& plane);
#endif
