//
//  Description of a right circular cylindrical, for use in KinKal intersection
//  original author: David Brown (LBN) 2023
//
#ifndef RecoGeom_Cylinder_hh
#define RecoGeom_Cylinder_hh
#include "Offline/RecoGeom/inc/Surface.hh"
namespace mu2e {
  namespace RecoGeom {
    class Cylinder : public Surface {
      public:
        // construct from necessary parameters
        Cylinder(XYZVectorD const& axis, XYZVectorD const& center, double radius, double halflen ) : axis_(axis), center_(center), radius_(radius), halflen_(halflen), radius2_(radius*radius) {}
        // surface interface

        bool onSurface(XYZVectorD const& point, double tol=1e-8) const override;
        IntersectFlag intersect(Ray const& ray,double& dist, double tol=1e-8) const override;
        XYZVectorD normal(XYZVectorD const& point) const override;
        // cylinder-specific interface
        auto const& axis() const { return axis_; }
        auto const& center() const { return center_; }
        double radius() const { return radius_; }
        double halfLength() const { return halflen_; }
      private:
        XYZVectorD axis_; // symmetry axis of the cylinder
        XYZVectorD center_; // geometric center
        double radius_; // transverse radius
        double halflen_; // half length
        double radius2_; // squared radius (cache);
    };
  }
}
std::ostream& operator <<(std::ostream& ost, mu2e::RecoGeom::Cylinder const& cyl);
#endif
