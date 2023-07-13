#include "Offline/RecoGeom/inc/Rectangle.hh"
#include "Math/VectorUtil.h"
using namespace ROOT::Math::VectorUtil;
namespace mu2e {
  namespace RecoGeom {
    Rectangle::Rectangle(XYZVectorD const& norm, XYZVectorD const& center, XYZVectorD const& uaxis, double uhalflen, double vhalflen) :
      Plane(norm,center) , uhalflen_(uhalflen), vhalflen_(vhalflen), udir_(uaxis.Unit()){
          // check that U is perpendicular
          if(udir_.Dot(normal()) > 1e-10) throw std::invalid_argument("U direction not perpendicular to normal");
          // V direction is implicit
          vdir_ = normal().Cross(udir_);
        }


    bool Rectangle::inBounds(XYZVectorD const& point, double tol) const {
      auto rvec = point - center();
      double udist = rvec.Dot(udir_);
      double vdist = rvec.Dot(vdir_);
      return onRectangle(udist,vdist);
    }
  }
}

std::ostream& operator <<(std::ostream& ost, mu2e::RecoGeom::Rectangle const& rect) {
  mu2e::RecoGeom::Plane const& plane = static_cast<mu2e::RecoGeom::Plane const&>(rect);
  ost << "Rectangle in " << plane << " U direction " << rect.uDirection() << " V direction " << rect.vDirection()
    << " U, V half-lengths " << rect.uHalfLength() << " , " << rect.vHalfLength();
  return ost;
}
