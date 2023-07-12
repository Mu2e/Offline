//
//  Description of an annular planar section
//  original author: David Brown (LBN) 2023
//
#ifndef RecoGeom_Annulus_hh
#define RecoGeom_Annulus_hh
#include "Offline/RecoGeom/inc/Plane.hh"
namespace mu2e {
  namespace RecoGeom {
    class Annulus : public Plane {
      public:
        // construct from necessary parameters
        Annulus(XYZVectorD const& norm, XYZVectorD const& center, double innerrad, double outerrad) : Plane(norm,center) {
          radii_[0] = innerrad;
          radii_[1] = outerrad;
        }
        auto innerRadius() const { return radii_[0]; }
        auto outerRadius() const { return radii_[1]; }
        bool onAnnulus(double radius) const { return radius >= radii_[0] && radius <= radii_[1]; }
        bool onSurface(XYZVectorD const& point, double tol=1e-8) const override;
        IntersectFlag intersect(Ray const& ray,double& dist, double tol=1e-8) const override;
      private:
        double radii_[2]; // inner and outer radii
    };
  }
}
std::ostream& operator <<(std::ostream& ost, mu2e::RecoGeom::Annulus const& annulus);

#endif
